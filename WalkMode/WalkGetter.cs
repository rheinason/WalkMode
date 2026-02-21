using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using Rhino;
using Rhino.Display;
using Rhino.DocObjects;
using Rhino.DocObjects.Tables;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using Rhino.Input;
using Rhino.Input.Custom;

namespace WalkMode
{
    internal sealed class WalkGetter : GetPoint
    {
        // ── P/Invoke ──────────────────────────────────────────────────────────
        [DllImport("user32.dll")]  static extern bool   SetCursorPos(int x, int y);
        // Win32 ShowCursor is reference-counted; Cursor.Hide() only calls it once
        // which may not be enough if Rhino has incremented the count internally.
        [DllImport("user32.dll", EntryPoint = "ShowCursor")]
        static extern int ShowCursorWin32(bool bShow);
        [DllImport("user32.dll")]  static extern bool   ClientToScreen(IntPtr hWnd, ref System.Drawing.Point lpPoint);
        [DllImport("user32.dll")]  static extern IntPtr SetWindowsHookEx(int idHook, LowLevelProc lpfn, IntPtr hMod, uint dwThreadId);
        [DllImport("user32.dll")]  static extern bool   UnhookWindowsHookEx(IntPtr hhk);
        [DllImport("user32.dll")]  static extern IntPtr CallNextHookEx(IntPtr hhk, int nCode, IntPtr wParam, IntPtr lParam);
        [DllImport("kernel32.dll")] static extern IntPtr GetModuleHandle(string lpModuleName);
        // Multimedia timer resolution — lets WinForms timer fire at 1 ms granularity
        [DllImport("winmm.dll")]   static extern int    timeBeginPeriod(int uPeriod);
        [DllImport("winmm.dll")]   static extern int    timeEndPeriod(int uPeriod);

        // Low-level hook IDs and message codes
        const int WH_KEYBOARD_LL = 13;
        const int WH_MOUSE_LL    = 14;
        const int WM_KEYDOWN     = 0x0100;
        const int WM_KEYUP       = 0x0101;
        const int WM_MOUSEMOVE   = 0x0200;
        const int WM_MOUSEWHEEL  = 0x020A;

        delegate IntPtr LowLevelProc(int nCode, IntPtr wParam, IntPtr lParam);

        [StructLayout(LayoutKind.Sequential)]
        struct POINT { public int x, y; }

        [StructLayout(LayoutKind.Sequential)]
        struct MSLLHOOKSTRUCT
        {
            public POINT  pt;
            public int    mouseData;
            public int    flags, time;
            public UIntPtr dwExtraInfo;
        }

        // ── Context ───────────────────────────────────────────────────────────
        readonly RhinoDoc      _doc;
        readonly RhinoView     _view;
        readonly RhinoViewport _vp;

        // ── Camera state ──────────────────────────────────────────────────────
        double  _yaw, _pitch;
        Point3d _camPos;

        // ── Cancel restore ────────────────────────────────────────────────────
        readonly ViewportInfo _savedViewport;

        // ── Unit scale (doc units per metre) ─────────────────────────────────
        double _unitsPerMeter;

        // ── Speed (doc units / s) ─────────────────────────────────────────────
        double _baseSpeed;

        // ── Cursor / screen centre ────────────────────────────────────────────
        System.Drawing.Point     _screenCenter;
        System.Drawing.Point     _lastMousePt;   // last known cursor pos; hook measures delta from this
        System.Drawing.Rectangle _lastBounds;
        // No cursor-hide count needed — we force the Win32 display counter
        // to exactly -1 every tick (hidden) and back to 0 in Cleanup (visible).

        // ── Frame timing ──────────────────────────────────────────────────────
        readonly Stopwatch _sw = Stopwatch.StartNew();
        long _lastTickMs;

        // ── Input state — written by hook callbacks, read by physics tick ─────
        // Accumulated mouse deltas (hook writes, timer reads+clears)
        int _mouseDx, _mouseDy;
        // Held-down keys
        bool _kW, _kA, _kS, _kD, _kQ, _kE, _kShift, _kAlt;
        // One-shot edges (set on keydown in hook, cleared after being consumed)
        bool _tabEdge, _vEdge, _spaceEdge;
        // Previous-down state for edge detection inside the hook
        bool _kTabDown, _kVDown, _kSpaceDown;

        // ── Gravity ───────────────────────────────────────────────────────────
        bool   _gravityEnabled;
        double _verticalVelocity;
        bool   _grounded;
        // Floor probe decimation — raycasts run at ~30 Hz, not every tick
        double _gravityProbeAccum;
        double _lastFloorZ;
        bool   _lastFloorValid;

        // ── Teleport ──────────────────────────────────────────────────────────
        bool     _teleporting;
        Point3d  _teleportOrigin, _teleportTarget;
        DateTime _teleportStart;
        double   _teleportDuration;   // distance-based, set in BeginTeleport

        // ── Scene geometry ────────────────────────────────────────────────────
        // All geometry is stored as meshes (render meshes extracted from Breps).
        // MeshRay needs no BVH warmup, so startup and first teleport are instant.
        Mesh[]   _meshGeometry;
        bool     _geometryDirty;
        DateTime _lastGeometryRebuild; // debounce: rebuild at most every 500 ms

        // ── Render throttle ──────────────────────────────────────────────────
        bool _redrawPending;
        long _lastRedrawDoneMs;

        // ── Near clip plane (doc units) — prevents walk-through clipping ────
        double _nearClip;

        // ── Doc event delegates (stored for unsubscribe) ──────────────────────
        EventHandler<RhinoObjectEventArgs>                 _objHandler;
        EventHandler<RhinoReplaceObjectEventArgs>          _replaceHandler;
        EventHandler<RhinoModifyObjectAttributesEventArgs> _attrHandler;
        EventHandler<LayerTableEventArgs>                  _layerHandler;

        // ── Low-level hooks ───────────────────────────────────────────────────
        LowLevelProc _hookProc;   // MUST be a field — prevents GC while the hook is live
        IntPtr       _kHook = IntPtr.Zero;
        IntPtr       _mHook = IntPtr.Zero;

        // ── Physics/render timer ──────────────────────────────────────────────
        System.Windows.Forms.Timer _timer;

        // ── Disposal guard ────────────────────────────────────────────────────
        bool _cleaned;

        // ── Sticky speed — persists across command invocations ────────────────
        static double s_lastSpeed = 0;

        // ─────────────────────────────────────────────────────────────────────
        public WalkGetter(RhinoDoc doc, RhinoView view)
        {
            _doc  = doc;
            _view = view;
            _vp   = view.ActiveViewport;

            _savedViewport = new ViewportInfo(_vp);

            _unitsPerMeter = RhinoMath.UnitScale(UnitSystem.Meters, _doc.ModelUnitSystem);

            // Speed: reuse last session value if available, otherwise derive from
            // camera distance (scene-aware) and clamp to a physical range.
            if (s_lastSpeed > 0)
            {
                _baseSpeed = s_lastSpeed;
            }
            else
            {
                double camDist = _vp.CameraLocation.DistanceTo(_vp.CameraTarget);
                _baseSpeed = Math.Max( 1.5 * _unitsPerMeter,
                             Math.Min(camDist / 50.0,
                                      50.0 * _unitsPerMeter));
            }

            var dir = _vp.CameraDirection;
            dir.Unitize();
            _yaw   = Math.Atan2(dir.X, dir.Y);
            _pitch = Math.Asin(Math.Max(-1.0, Math.Min(1.0, dir.Z)));
            _camPos = _vp.CameraLocation;
            AcceptNothing(true);
            PermitObjectSnap(false);
            UpdatePrompt();

            BuildGeometryCache();
            _objHandler     = (s, e) => _geometryDirty = true;
            _replaceHandler = (s, e) => _geometryDirty = true;
            _attrHandler    = (s, e) => _geometryDirty = true;
            _layerHandler   = (s, e) => _geometryDirty = true;
            RhinoDoc.AddRhinoObject         += _objHandler;
            RhinoDoc.DeleteRhinoObject      += _objHandler;
            RhinoDoc.ReplaceRhinoObject     += _replaceHandler;
            RhinoDoc.ModifyObjectAttributes += _attrHandler;
            RhinoDoc.LayerTableEvent        += _layerHandler;

            _lastBounds = _vp.Bounds;
            // Near clip: 10 cm equivalent — close enough to walk past walls
            // without clipping, far enough to avoid z-fighting.
            _nearClip = 0.1 * _unitsPerMeter;
            ComputeScreenCenter();

            // Install global low-level keyboard + mouse hooks.
            // These intercept input at the OS level, before any window (including
            // Rhino's command line) can see it.  Returning IntPtr(-1) consumes the
            // event so it never types "W" into the command line.
            // The delegate MUST be stored in a field to prevent the GC from
            // collecting it while the hook is installed.
            _hookProc = HookCallback;
            using (var proc = Process.GetCurrentProcess())
            using (var mod  = proc.MainModule)
            {
                var hMod = GetModuleHandle(mod.ModuleName);
                _kHook = SetWindowsHookEx(WH_KEYBOARD_LL, _hookProc, hMod, 0);
                _mHook = SetWindowsHookEx(WH_MOUSE_LL,    _hookProc, hMod, 0);
            }

            _lastTickMs = _sw.ElapsedMilliseconds;

            // 4 ms (250 Hz) for low-latency input.  Previously this caused lag
            // because _view.Redraw() blocked the message pump.  Now we use
            // InvalidateRect which returns instantly, so ticks are ~0.1 ms each.
            // timeBeginPeriod(1) ensures the OS actually fires at 4 ms granularity.
            timeBeginPeriod(1);
            _timer = new System.Windows.Forms.Timer { Interval = 4 };
            _timer.Tick += OnTimerTick;
            _timer.Start();
        }

        // ── Geometry cache ────────────────────────────────────────────────────
        void BuildGeometryCache()
        {
            var meshes = new List<Mesh>();
            CollectGeometry(_doc.Objects, Transform.Identity, meshes);
            _meshGeometry  = meshes.ToArray();
            _geometryDirty = false;
        }

        // Recursively collects geometry as meshes, expanding block instances.
        // Breps/Surfaces contribute their pre-built render meshes (RhinoObject.GetMeshes)
        // — no BVH build, no warmup needed.  Falls back to fast meshing for block
        // members whose render meshes haven't been built yet.
        void CollectGeometry(IEnumerable<RhinoObject> objects, Transform xform,
                             List<Mesh> meshes)
        {
            foreach (var obj in objects)
            {
                if (!obj.Visible) continue;

                if (obj is InstanceObject inst)
                {
                    var childXform = xform * inst.InstanceXform;
                    CollectGeometry(inst.InstanceDefinition.GetObjects(), childXform, meshes);
                    continue;
                }

                var geom = obj.Geometry;
                Mesh[] source = null;

                if (geom is Mesh m)
                {
                    source = new[] { m };
                }
                else if (geom is Brep || geom is Surface)
                {
                    source = obj.GetMeshes(MeshType.Render);
                    // Fallback for block-definition objects that haven't been rendered yet.
                    if ((source == null || source.Length == 0) && geom is Brep brep)
                        source = Mesh.CreateFromBrep(brep, MeshingParameters.FastRenderMesh);
                }

                if (source == null) continue;
                foreach (var rm in source)
                {
                    if (rm == null) continue;
                    if (xform.IsIdentity)
                    {
                        meshes.Add(rm);
                    }
                    else
                    {
                        var dup = (Mesh)rm.Duplicate();
                        dup.Transform(xform);
                        meshes.Add(dup);
                    }
                }
            }
        }

        // Closest mesh hit; null = no hit.
        Point3d? ShootRay(Ray3d ray)
        {
            // Debounce: rebuild at most every 500 ms so a burst of doc-change events
            // (layer toggle, object edit) doesn't cause a synchronous stall mid-walk.
            if (_geometryDirty)
            {
                var now = DateTime.UtcNow;
                if ((now - _lastGeometryRebuild).TotalSeconds >= 0.5)
                {
                    BuildGeometryCache();
                    _lastGeometryRebuild = now;
                }
            }

            Point3d? closest     = null;
            double   closestDist = double.MaxValue;

            if (_meshGeometry != null)
            {
                foreach (var mesh in _meshGeometry)
                {
                    double t = Intersection.MeshRay(mesh, ray);
                    if (t >= 0 && t < closestDist)
                    {
                        closestDist = t;
                        closest = ray.Position + ray.Direction * t;
                    }
                }
            }

            return closest;
        }

        void ComputeScreenCenter()
        {
            // Use _vp.Bounds (viewport rect in view-window client coordinates),
            // not _view.ClientRectangle, so the center lands in the active
            // viewport even when the window contains multiple viewports.
            var bounds = _vp.Bounds;
            var pt = new System.Drawing.Point(
                bounds.X + bounds.Width  / 2,
                bounds.Y + bounds.Height / 2);
            if (_view.Handle != IntPtr.Zero)
                ClientToScreen(_view.Handle, ref pt);
            _screenCenter = pt;
            _lastMousePt  = pt;   // reset reference to avoid a delta jump on resize
        }

        void OnWheelDelta(int delta)
        {
            if (delta > 0) _baseSpeed *= 1.2;
            else           _baseSpeed /= 1.2;
            _baseSpeed = Math.Max(0.1 * _unitsPerMeter, Math.Min(_baseSpeed, 100.0 * _unitsPerMeter));
            s_lastSpeed = _baseSpeed;
            UpdatePrompt();
        }

        // ── Low-level hook callback ───────────────────────────────────────────
        // Called on every keyboard / mouse event system-wide, on the UI thread
        // (because the hook was installed from the UI thread and the message pump
        // is running inside GetPoint.Get()).
        // Returning IntPtr(-1) consumes the event — Rhino never sees it.
        IntPtr HookCallback(int nCode, IntPtr wParam, IntPtr lParam)
        {
            if (nCode < 0)
                return CallNextHookEx(IntPtr.Zero, nCode, wParam, lParam);

            int msg = (int)wParam;

            // ── Keyboard ──────────────────────────────────────────────────────
            if (msg == WM_KEYDOWN || msg == WM_KEYUP)
            {
                bool down = (msg == WM_KEYDOWN);
                var  key  = (Keys)Marshal.ReadInt32(lParam);
                bool consume = true;

                switch (key)
                {
                    case Keys.W:                               _kW     = down; break;
                    case Keys.A:                               _kA     = down; break;
                    case Keys.S:                               _kS     = down; break;
                    case Keys.D:                               _kD     = down; break;
                    case Keys.Q:                               _kQ     = down; break;
                    case Keys.E:                               _kE     = down; break;
                    case Keys.LShiftKey: case Keys.RShiftKey:  _kShift = down; break;
                    case Keys.LMenu:     case Keys.RMenu:      _kAlt   = down; break;
                    case Keys.Tab:
                        if (down && !_kTabDown) _tabEdge = true;
                        _kTabDown = down;
                        break;
                    case Keys.V:
                        if (down && !_kVDown) _vEdge = true;
                        _kVDown = down;
                        break;
                    case Keys.Space:
                        if (down && !_kSpaceDown) _spaceEdge = true;
                        _kSpaceDown = down;
                        break;
                    default:
                        consume = false;    // pass Esc, Enter, etc. through to Rhino
                        break;
                }

                if (consume) return new IntPtr(-1);
            }
            // ── Mouse move ────────────────────────────────────────────────────
            else if (msg == WM_MOUSEMOVE)
            {
                var data = Marshal.PtrToStructure<MSLLHOOKSTRUCT>(lParam);
                int dx = data.pt.x - _lastMousePt.X;
                int dy = data.pt.y - _lastMousePt.Y;
                _lastMousePt = new System.Drawing.Point(data.pt.x, data.pt.y);
                if (dx != 0 || dy != 0)
                {
                    // Just accumulate raw deltas — no Rhino API calls here.
                    // The timer tick consumes these at 250 Hz (4 ms), applies
                    // the camera once, and invalidates once.  Keeping the hook
                    // ultra-light prevents it from starving WM_TIMER messages.
                    _mouseDx += dx;
                    _mouseDy += dy;

                    // Warp cursor back to center so it never hits the screen edge.
                    // Echo WM_MOUSEMOVE has dx=dy=0 → skipped, no cascade.
                    SetCursorPos(_screenCenter.X, _screenCenter.Y);
                    _lastMousePt = _screenCenter;
                }
                return new IntPtr(-1);   // consume — Rhino doesn't see cursor move
            }
            // ── Scroll wheel ──────────────────────────────────────────────────
            else if (msg == WM_MOUSEWHEEL)
            {
                var data  = Marshal.PtrToStructure<MSLLHOOKSTRUCT>(lParam);
                int delta = (short)(data.mouseData >> 16);
                OnWheelDelta(delta);
                return new IntPtr(-1);
            }

            return CallNextHookEx(IntPtr.Zero, nCode, wParam, lParam);
        }

        // ── Physics timer tick (250 Hz) ──────────────────────────────────────
        // All input (mouse deltas, keyboard, gravity) is consumed here in one
        // place.  ApplyCamera runs every tick something changed (~0.01 ms);
        // Redraw is throttled so the synchronous render can't back up the queue.
        void OnTimerTick(object sender, EventArgs e)
        {
            // Force the Win32 cursor-display counter to exactly -1 (hidden).
            {
                int c = ShowCursorWin32(false);
                while (c > -1) c = ShowCursorWin32(false);
                while (c < -1) c = ShowCursorWin32(true);
            }

            long   nowMs = _sw.ElapsedMilliseconds;
            double dt    = Math.Min((nowMs - _lastTickMs) / 1000.0, 0.1);
            _lastTickMs  = nowMs;

            // Refresh screen centre if the viewport was resized or moved
            var currentBounds = _vp.Bounds;
            if (currentBounds != _lastBounds)
            {
                _lastBounds = currentBounds;
                ComputeScreenCenter();
            }

            // ── Consume accumulated mouse deltas ─────────────────────────────
            int mdx = _mouseDx, mdy = _mouseDy;
            _mouseDx = 0; _mouseDy = 0;
            if (mdx != 0 || mdy != 0)
            {
                const double sens = 0.002;
                _yaw   += mdx * sens;
                _pitch -= mdy * sens;
                const double maxPitch = 89.0 * Math.PI / 180.0;
                if (_pitch >  maxPitch) _pitch =  maxPitch;
                if (_pitch < -maxPitch) _pitch = -maxPitch;
                _redrawPending = true;
            }

            bool moved = ProcessKeyboardMovement(dt);
            if (_gravityEnabled && !_teleporting) moved |= ProcessGravity(dt);
            if (_teleporting)                     { ProcessTeleport(); moved = true; }
            ProcessToggles();

            if (moved) _redrawPending = true;

            // Apply camera state every tick something changed (cheap).
            // Redraw only when enough time has passed since the last frame
            // completed — this leaves the message pump free to dispatch
            // mouse hooks between frames while keeping a steady framerate.
            if (_redrawPending)
            {
                ApplyCamera();

                if (nowMs - _lastRedrawDoneMs >= 6)
                {
                    _redrawPending = false;
                    _view.Redraw();
                    _lastRedrawDoneMs = _sw.ElapsedMilliseconds;
                }
            }
        }

        // ── HUD draw ──────────────────────────────────────────────────────────
        protected override void OnDynamicDraw(GetPointDrawEventArgs e)
        {
            base.OnDynamicDraw(e);
            DrawCrosshair(e);
            DrawStatus(e);
        }

        void DrawStatus(GetPointDrawEventArgs e)
        {
            double speedMs = _baseSpeed / _unitsPerMeter;
            string grav    = _gravityEnabled ? "Gravity ON" : "Gravity OFF";
            string text    = $"{grav}   Speed {speedMs:F1} m/s";
            int    y       = _vp.Bounds.Height - 24;
            e.Display.Draw2dText(text, System.Drawing.Color.White,
                new Point2d(10, y), false, 13);
        }

        // ── Movement / physics ────────────────────────────────────────────────
        bool ProcessKeyboardMovement(double dt)
        {
            double speed = _baseSpeed;
            if (_kShift) speed *= 5.0;
            if (_kAlt)   speed *= 0.2;

            var fwd   = new Vector3d( Math.Sin(_yaw),  Math.Cos(_yaw), 0);
            var right = new Vector3d( Math.Cos(_yaw), -Math.Sin(_yaw), 0);

            var move = Vector3d.Zero;
            if (_kW) move += fwd;
            if (_kS) move -= fwd;
            if (_kD) move += right;
            if (_kA) move -= right;

            if (!_gravityEnabled)
            {
                if (_kQ) move -= Vector3d.ZAxis;
                if (_kE) move += Vector3d.ZAxis;
            }

            if (move.SquareLength > 0)
            {
                move.Unitize();
                _camPos = _camPos + move * (speed * dt);
                return true;
            }
            return false;
        }

        bool ProcessGravity(double dt)
        {
            double g         = -9.8  * _unitsPerMeter;
            double terminal  = -50.0 * _unitsPerMeter;
            double eyeHeight =  1.75 * _unitsPerMeter;

            double prevZ = _camPos.Z;

            _verticalVelocity += g * dt;
            if (_verticalVelocity < terminal) _verticalVelocity = terminal;
            _camPos.Z += _verticalVelocity * dt;

            // Probe the floor at ~30 Hz rather than every tick.
            // Between probes, snap against the cached floor height.
            const double ProbeHz = 30.0;
            _gravityProbeAccum += dt;
            if (_gravityProbeAccum >= 1.0 / ProbeHz)
            {
                _gravityProbeAccum = 0;
                var floor = ShootRay(new Ray3d(_camPos, -Vector3d.ZAxis));
                _lastFloorValid = floor.HasValue;
                if (floor.HasValue) _lastFloorZ = floor.Value.Z;
            }

            _grounded = false;
            if (_lastFloorValid)
            {
                double minZ = _lastFloorZ + eyeHeight;
                if (_camPos.Z <= minZ)
                {
                    _camPos.Z         = minZ;
                    _verticalVelocity = 0;
                    _grounded         = true;
                }
            }

            // Only dirty the camera when Z actually moved — avoids a 60 Hz
            // full-viewport redraw while the player is standing still on a floor.
            return Math.Abs(_camPos.Z - prevZ) > 1e-9;
        }

        void ProcessTeleport()
        {
            double elapsed = (DateTime.UtcNow - _teleportStart).TotalSeconds;
            double t       = Math.Min(elapsed / _teleportDuration, 1.0);
            double smooth  = t * t * (3.0 - 2.0 * t);

            _camPos = _teleportOrigin + (_teleportTarget - _teleportOrigin) * smooth;

            if (t >= 1.0)
            {
                _teleporting      = false;
                _verticalVelocity = 0;
            }
        }

        void BeginTeleport()
        {
            // Use the GPU depth buffer — reads whatever Rhino has already rendered,
            // so blocks, instances, meshes, NURBS, everything works automatically.
            // No geometry cache needed, no ray-cast stall.
            Point3d hit;
            using (var zBuf = new ZBufferCapture(_vp))
            {
                int cx = _vp.Bounds.Width  / 2;
                int cy = _vp.Bounds.Height / 2;
                hit = zBuf.WorldPointAt(cx, cy);
            }

            if (!hit.IsValid) return;
            double dist = _camPos.DistanceTo(hit);
            if (dist < 0.1 * _unitsPerMeter) return;   // pointing at sky or self

            _teleportOrigin   = _camPos;
            _teleportTarget   = new Point3d(hit.X, hit.Y, hit.Z + 1.75 * _unitsPerMeter);
            _teleportStart    = DateTime.UtcNow;
            // Scale duration with distance: ~60 ms for close hops, up to 150 ms.
            _teleportDuration = Math.Max(0.06, Math.Min(0.15, dist / (20.0 * _unitsPerMeter)));
            _teleporting      = true;
        }

        void ProcessToggles()
        {
            // Tab → toggle gravity
            if (_tabEdge)
            {
                _tabEdge = false;
                _gravityEnabled = !_gravityEnabled;
                if (!_gravityEnabled)
                {
                    _verticalVelocity = 0;
                    _grounded = false;
                }
                UpdatePrompt();
            }

            // V → jump (only when confirmed grounded by floor probe)
            if (_vEdge)
            {
                _vEdge = false;
                if (_gravityEnabled && _grounded)
                {
                    _verticalVelocity = 5.0 * _unitsPerMeter;
                    _grounded = false;
                }
            }

            // Space → teleport
            if (_spaceEdge)
            {
                _spaceEdge = false;
                if (!_teleporting) BeginTeleport();
            }
        }

        void UpdatePrompt()
        {
            double speedMs = _baseSpeed / _unitsPerMeter;
            string grav    = _gravityEnabled ? "Gravity ON" : "Gravity OFF";
            SetCommandPrompt(
                $"Walk  WASD·move  Mouse·look  Tab·{grav}  V·jump  Space·teleport  Speed·{speedMs:F1}m/s·(scroll)  Esc·exit");
        }

        Vector3d ReconstructCameraDirection()
        {
            return new Vector3d(
                Math.Cos(_pitch) * Math.Sin(_yaw),
                Math.Cos(_pitch) * Math.Cos(_yaw),
                Math.Sin(_pitch));
        }

        void ApplyCamera()
        {
            var dir = ReconstructCameraDirection();
            _vp.SetCameraLocations(_camPos + dir, _camPos);

            // Force a walkthrough-friendly near clip plane.  Rhino's auto-
            // calculation bases near/far on scene extents, which in mm or cm
            // files pushes the near plane so far out that nearby walls clip.
            var vi = new ViewportInfo(_vp);
            if (vi.FrustumNear > _nearClip * 2)
            {
                vi.SetFrustumNearFar(_nearClip, vi.FrustumFar);
                _vp.SetViewProjection(vi, false);
            }
        }

        // ── Crosshair HUD ─────────────────────────────────────────────────────
        void DrawCrosshair(GetPointDrawEventArgs e)
        {
            int cx = _vp.Bounds.Width  / 2;
            int cy = _vp.Bounds.Height / 2;

            const int   arm = 12;
            const int   gap = 4;
            const float thk = 1.5f;
            var white = System.Drawing.Color.White;

            e.Display.Draw2dLine(
                new System.Drawing.PointF(cx,       cy - gap),
                new System.Drawing.PointF(cx,       cy - gap - arm), white, thk);
            e.Display.Draw2dLine(
                new System.Drawing.PointF(cx,       cy + gap),
                new System.Drawing.PointF(cx,       cy + gap + arm), white, thk);
            e.Display.Draw2dLine(
                new System.Drawing.PointF(cx - gap, cy),
                new System.Drawing.PointF(cx - gap - arm, cy),       white, thk);
            e.Display.Draw2dLine(
                new System.Drawing.PointF(cx + gap, cy),
                new System.Drawing.PointF(cx + gap + arm, cy),       white, thk);
            e.Display.Draw2dLine(
                new System.Drawing.PointF(cx - 1, cy),
                new System.Drawing.PointF(cx + 1, cy), white, 2.0f);
        }

        // ── Entry point ───────────────────────────────────────────────────────
        internal GetResult RunWalk()
        {
            if (!_vp.IsPerspectiveProjection)
            {
                RhinoApp.WriteLine("WalkMode requires a perspective viewport.");
                return GetResult.Cancel;
            }

            // Position cursor at viewport center before Get() starts.
            // Cursor hiding is deferred to the first timer tick — Get() internally
            // calls ShowCursor(true) which undoes an early Cursor.Hide(), so we
            // drain the Win32 reference count from inside the running message pump.
            SetCursorPos(_screenCenter.X, _screenCenter.Y);

            GetResult result;
            try
            {
                result = Get();
            }
            finally
            {
                Cleanup();
            }

            if (result == GetResult.Cancel)
                RestoreCamera();

            return result;
        }

        void RestoreCamera()
        {
            _vp.SetViewProjection(_savedViewport, false);
        }

        // ── Cleanup ───────────────────────────────────────────────────────────
        void Cleanup()
        {
            if (_cleaned) return;
            _cleaned = true;

            // Persist speed so the next run starts where this one left off.
            s_lastSpeed = _baseSpeed;

            // Restore default OS timer resolution, then stop the physics timer.
            timeEndPeriod(1);
            if (_timer != null) { _timer.Stop(); _timer.Dispose(); _timer = null; }

            // Uninstall global hooks — stop the callback from firing during teardown.
            if (_kHook != IntPtr.Zero) { UnhookWindowsHookEx(_kHook); _kHook = IntPtr.Zero; }
            if (_mHook != IntPtr.Zero) { UnhookWindowsHookEx(_mHook); _mHook = IntPtr.Zero; }

            // Unsubscribe doc events
            RhinoDoc.AddRhinoObject         -= _objHandler;
            RhinoDoc.DeleteRhinoObject      -= _objHandler;
            RhinoDoc.ReplaceRhinoObject     -= _replaceHandler;
            RhinoDoc.ModifyObjectAttributes -= _attrHandler;
            RhinoDoc.LayerTableEvent        -= _layerHandler;

            // Force cursor visible: restore Win32 display counter to exactly 0.
            {
                int c = ShowCursorWin32(true);
                while (c < 0) c = ShowCursorWin32(true);
                while (c > 0) c = ShowCursorWin32(false);
            }

            _meshGeometry = null;
        }

        protected override void Dispose(bool disposing)
        {
            if (disposing) Cleanup();
            base.Dispose(disposing);
        }
    }
}
