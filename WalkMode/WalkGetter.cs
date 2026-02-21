using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using Rhino;
using Rhino.Commands;
using Rhino.Display;
using Rhino.DocObjects;
using Rhino.DocObjects.Tables;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace WalkMode
{
    internal sealed class WalkGetter : IDisposable
    {
        // ── P/Invoke ──────────────────────────────────────────────────────────
        [DllImport("user32.dll")]  static extern bool   SetCursorPos(int x, int y);
        [DllImport("user32.dll", EntryPoint = "ShowCursor")]
        static extern int ShowCursorWin32(bool bShow);
        [DllImport("user32.dll")]  static extern bool   ClientToScreen(IntPtr hWnd, ref System.Drawing.Point lpPoint);
        [DllImport("user32.dll")]  static extern IntPtr SetWindowsHookEx(int idHook, LowLevelProc lpfn, IntPtr hMod, uint dwThreadId);
        [DllImport("user32.dll")]  static extern bool   UnhookWindowsHookEx(IntPtr hhk);
        [DllImport("user32.dll")]  static extern IntPtr CallNextHookEx(IntPtr hhk, int nCode, IntPtr wParam, IntPtr lParam);
        [DllImport("kernel32.dll")] static extern IntPtr GetModuleHandle(string lpModuleName);
        [DllImport("winmm.dll")]   static extern int    timeBeginPeriod(int uPeriod);
        [DllImport("winmm.dll")]   static extern int    timeEndPeriod(int uPeriod);

        // Message pump P/Invokes
        [DllImport("user32.dll")]
        static extern bool PeekMessage(out MSG lpMsg, IntPtr hWnd, uint wMsgFilterMin, uint wMsgFilterMax, uint wRemoveMsg);
        [DllImport("user32.dll")]
        static extern bool TranslateMessage(ref MSG lpMsg);
        [DllImport("user32.dll")]
        static extern IntPtr DispatchMessage(ref MSG lpMsg);
        [DllImport("user32.dll")]
        static extern bool WaitMessage();

        const uint PM_REMOVE = 0x0001;

        [StructLayout(LayoutKind.Sequential)]
        struct MSG
        {
            public IntPtr hwnd;
            public uint   message;
            public IntPtr wParam;
            public IntPtr lParam;
            public uint   time;
            public POINT  pt;
        }

        // Low-level hook IDs and message codes
        const int WH_KEYBOARD_LL = 13;
        const int WH_MOUSE_LL    = 14;
        const int WM_KEYDOWN     = 0x0100;
        const int WM_KEYUP       = 0x0101;
        const int WM_MOUSEMOVE   = 0x0200;
        const int WM_RBUTTONDOWN = 0x0204;
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
        System.Drawing.Point     _lastMousePt;
        System.Drawing.Rectangle _lastBounds;

        // ── Frame timing ──────────────────────────────────────────────────────
        readonly Stopwatch _sw = Stopwatch.StartNew();
        long _lastTickMs;

        // ── Input state — written by hook callbacks, read by physics tick ─────
        int _mouseDx, _mouseDy;
        bool _kW, _kA, _kS, _kD, _kQ, _kE, _kShift, _kAlt;
        bool _tabEdge, _vEdge, _spaceEdge;
        bool _kTabDown, _kVDown, _kSpaceDown;

        // ── Exit flags (replace GetPoint's Get() return) ─────────────────────
        bool _exitRequested;
        bool _cancelled;

        // ── Gravity ───────────────────────────────────────────────────────────
        bool   _gravityEnabled;
        double _verticalVelocity;
        bool   _grounded;
        double _gravityProbeAccum;
        double _lastFloorZ;
        bool   _lastFloorValid;

        // ── Teleport ──────────────────────────────────────────────────────────
        bool     _teleporting;
        Point3d  _teleportOrigin, _teleportTarget;
        DateTime _teleportStart;
        double   _teleportDuration;

        // ── Scene geometry ────────────────────────────────────────────────────
        Mesh[]   _meshGeometry;
        bool     _geometryDirty;
        DateTime _lastGeometryRebuild;

        // ── Render state ────────────────────────────────────────────────────
        bool _redrawPending;

        // ── Near clip plane (doc units) ───────────────────────────────────────
        double _nearClip;

        // ── Doc event delegates ─────────────────────────────────────────────
        EventHandler<RhinoObjectEventArgs>                 _objHandler;
        EventHandler<RhinoReplaceObjectEventArgs>          _replaceHandler;
        EventHandler<RhinoModifyObjectAttributesEventArgs> _attrHandler;
        EventHandler<LayerTableEventArgs>                  _layerHandler;

        // ── Low-level hooks ───────────────────────────────────────────────────
        LowLevelProc _hookProc;
        IntPtr       _kHook = IntPtr.Zero;
        IntPtr       _mHook = IntPtr.Zero;

        // ── Physics/render timer ──────────────────────────────────────────────
        System.Windows.Forms.Timer _timer;

        // ── Display conduit ───────────────────────────────────────────────────
        WalkConduit _conduit;

        // ── Disposal guard ────────────────────────────────────────────────────
        bool _cleaned;

        // ── Sticky speed ──────────────────────────────────────────────────────
        static double s_lastSpeed = 0;

        // ── Display conduit for crosshair + status HUD ────────────────────────
        sealed class WalkConduit : DisplayConduit
        {
            readonly WalkGetter _owner;
            readonly Guid       _vpId;

            public WalkConduit(WalkGetter owner)
            {
                _owner = owner;
                _vpId  = owner._vp.Id;
            }

            protected override void DrawForeground(DrawEventArgs e)
            {
                if (e.Viewport.Id != _vpId) return;

                DrawCrosshair(e);
                DrawStatus(e);
            }

            void DrawCrosshair(DrawEventArgs e)
            {
                int cx = _owner._vp.Bounds.Width  / 2;
                int cy = _owner._vp.Bounds.Height / 2;

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

            void DrawStatus(DrawEventArgs e)
            {
                double speedMs = _owner._baseSpeed / _owner._unitsPerMeter;
                string grav    = _owner._gravityEnabled ? "Gravity ON" : "Gravity OFF";
                string text    = $"{grav}   Speed {speedMs:F1} m/s";
                int    y       = _owner._vp.Bounds.Height - 24;
                e.Display.Draw2dText(text, System.Drawing.Color.White,
                    new Point2d(10, y), false, 13);
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        public WalkGetter(RhinoDoc doc, RhinoView view)
        {
            _doc  = doc;
            _view = view;
            _vp   = view.ActiveViewport;

            _savedViewport = new ViewportInfo(_vp);

            _unitsPerMeter = RhinoMath.UnitScale(UnitSystem.Meters, _doc.ModelUnitSystem);

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
            _nearClip = 0.1 * _unitsPerMeter;
            ComputeScreenCenter();

            _hookProc = HookCallback;
            using (var proc = Process.GetCurrentProcess())
            using (var mod  = proc.MainModule)
            {
                var hMod = GetModuleHandle(mod.ModuleName);
                _kHook = SetWindowsHookEx(WH_KEYBOARD_LL, _hookProc, hMod, 0);
                _mHook = SetWindowsHookEx(WH_MOUSE_LL,    _hookProc, hMod, 0);
            }

            _lastTickMs = _sw.ElapsedMilliseconds;

            timeBeginPeriod(1);
            _timer = new System.Windows.Forms.Timer { Interval = 4 };
            _timer.Tick += OnTimerTick;
            _timer.Start();

            // Enable display conduit for HUD drawing
            _conduit = new WalkConduit(this);
            _conduit.Enabled = true;
        }

        // ── Geometry cache ────────────────────────────────────────────────────
        void BuildGeometryCache()
        {
            var meshes = new List<Mesh>();
            CollectGeometry(_doc.Objects, Transform.Identity, meshes);
            _meshGeometry  = meshes.ToArray();
            _geometryDirty = false;
        }

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

        Point3d? ShootRay(Ray3d ray)
        {
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
            var bounds = _vp.Bounds;
            var pt = new System.Drawing.Point(
                bounds.X + bounds.Width  / 2,
                bounds.Y + bounds.Height / 2);
            if (_view.Handle != IntPtr.Zero)
                ClientToScreen(_view.Handle, ref pt);
            _screenCenter = pt;
            _lastMousePt  = pt;
        }

        void OnWheelDelta(int delta)
        {
            if (delta > 0) _baseSpeed *= 1.2;
            else           _baseSpeed /= 1.2;
            _baseSpeed = Math.Max(0.1 * _unitsPerMeter, Math.Min(_baseSpeed, 100.0 * _unitsPerMeter));
            s_lastSpeed = _baseSpeed;
            _redrawPending = true;
        }

        // ── Low-level hook callback ───────────────────────────────────────────
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
                    case Keys.Escape:
                        if (down)
                        {
                            _cancelled     = true;
                            _exitRequested = true;
                        }
                        break;
                    case Keys.Return:
                        if (down)
                        {
                            _exitRequested = true;
                        }
                        break;
                    default:
                        consume = false;
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
                    _mouseDx += dx;
                    _mouseDy += dy;
                    SetCursorPos(_screenCenter.X, _screenCenter.Y);
                    _lastMousePt = _screenCenter;
                }
                return new IntPtr(-1);
            }
            // ── Right-click → accept view ─────────────────────────────────────
            else if (msg == WM_RBUTTONDOWN)
            {
                _exitRequested = true;
                return new IntPtr(-1);
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

            if (_redrawPending)
            {
                _redrawPending = false;
                ApplyCamera();
                AdjustNearClip();
                _view.Redraw();
            }
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
            Point3d hit;
            using (var zBuf = new ZBufferCapture(_vp))
            {
                int cx = _vp.Bounds.Width  / 2;
                int cy = _vp.Bounds.Height / 2;
                hit = zBuf.WorldPointAt(cx, cy);
            }

            if (!hit.IsValid) return;
            double dist = _camPos.DistanceTo(hit);
            if (dist < 0.1 * _unitsPerMeter) return;

            _teleportOrigin   = _camPos;
            _teleportTarget   = new Point3d(hit.X, hit.Y, hit.Z + 1.75 * _unitsPerMeter);
            _teleportStart    = DateTime.UtcNow;
            _teleportDuration = Math.Max(0.06, Math.Min(0.15, dist / (20.0 * _unitsPerMeter)));
            _teleporting      = true;
        }

        void ProcessToggles()
        {
            if (_tabEdge)
            {
                _tabEdge = false;
                _gravityEnabled = !_gravityEnabled;
                if (!_gravityEnabled)
                {
                    _verticalVelocity = 0;
                    _grounded = false;
                }
                _redrawPending = true;
            }

            if (_vEdge)
            {
                _vEdge = false;
                if (_gravityEnabled && _grounded)
                {
                    _verticalVelocity = 5.0 * _unitsPerMeter;
                    _grounded = false;
                }
            }

            if (_spaceEdge)
            {
                _spaceEdge = false;
                if (!_teleporting) BeginTeleport();
            }
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
        }

        void AdjustNearClip()
        {
            var vi = new ViewportInfo(_vp);
            if (vi.FrustumNear > _nearClip * 2)
            {
                vi.SetFrustumNearFar(_nearClip, vi.FrustumFar);
                _vp.SetViewProjection(vi, false);
            }
        }

        // ── Entry point ───────────────────────────────────────────────────────
        internal Result RunWalk()
        {
            if (!_vp.IsPerspectiveProjection)
            {
                RhinoApp.WriteLine("WalkMode requires a perspective viewport.");
                return Result.Cancel;
            }

            SetCursorPos(_screenCenter.X, _screenCenter.Y);

            RhinoApp.SetCommandPrompt(
                "Walk  WASD·move  Mouse·look  Tab·gravity  V·jump  Space·teleport  Scroll·speed  Esc·exit");

            try
            {
                // Manual message pump — lightest possible loop.
                // WaitMessage blocks efficiently until the OS has a message to deliver,
                // so CPU usage is zero when idle.  The timer, hooks, and Rhino's own
                // internal messages all fire normally within this pump.
                while (!_exitRequested)
                {
                    while (PeekMessage(out var msg, IntPtr.Zero, 0, 0, PM_REMOVE))
                    {
                        TranslateMessage(ref msg);
                        DispatchMessage(ref msg);
                    }
                    WaitMessage();
                }
            }
            finally
            {
                Cleanup();
            }

            if (_cancelled)
            {
                RestoreCamera();
                return Result.Cancel;
            }

            return Result.Success;
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

            s_lastSpeed = _baseSpeed;

            // Disable display conduit
            if (_conduit != null) { _conduit.Enabled = false; _conduit = null; }

            timeEndPeriod(1);
            if (_timer != null) { _timer.Stop(); _timer.Dispose(); _timer = null; }

            if (_kHook != IntPtr.Zero) { UnhookWindowsHookEx(_kHook); _kHook = IntPtr.Zero; }
            if (_mHook != IntPtr.Zero) { UnhookWindowsHookEx(_mHook); _mHook = IntPtr.Zero; }

            RhinoDoc.AddRhinoObject         -= _objHandler;
            RhinoDoc.DeleteRhinoObject      -= _objHandler;
            RhinoDoc.ReplaceRhinoObject     -= _replaceHandler;
            RhinoDoc.ModifyObjectAttributes -= _attrHandler;
            RhinoDoc.LayerTableEvent        -= _layerHandler;

            // Force cursor visible
            {
                int c = ShowCursorWin32(true);
                while (c < 0) c = ShowCursorWin32(true);
                while (c > 0) c = ShowCursorWin32(false);
            }

            _meshGeometry = null;
        }

        public void Dispose()
        {
            Cleanup();
        }
    }
}
