# WalkMode

A fast first-person walkthrough camera plugin for Rhino 8. Uses a lightweight Win32 message pump instead of Rhino's GetPoint framework, giving smooth mouse-look and movement that matches native Rhino performance.

## Install

Search for **WalkMode** in Rhino's Package Manager (`_PackageManager`), or install via the Yak CLI:

```
"C:\Program Files\Rhino 8\System\Yak.exe" install WalkMode
```

## Usage

Run the `WalkMode` command in a perspective viewport.

### Controls

| Input | Action |
|---|---|
| **Mouse** | Look around |
| **W A S D** | Move forward / left / back / right |
| **Q / E** | Move down / up (fly mode only) |
| **Shift** | Sprint (5x speed) |
| **Alt** | Slow walk (0.2x speed) |
| **Scroll wheel** | Adjust base speed |
| **Tab** | Toggle gravity on/off |
| **V** | Jump (when gravity is on and grounded) |
| **Space** | Teleport to crosshair target |
| **Enter** or **Right-click** | Accept current view |
| **Escape** | Cancel and restore original view |

### Features

- **Gravity mode** with ground detection, jumping, and terminal velocity
- **Teleport** to any visible surface via depth buffer (works with all geometry types)
- **Adaptive near-clip** for mm/cm scenes — no walk-through clipping
- **Speed persistence** across command invocations
- **HUD** with crosshair, gravity status, and speed readout

## Building from source

Requires the .NET SDK and Rhino 8 installed.

```
dotnet build
```

Outputs both `net48` and `net7.0-windows` builds to `WalkMode/bin/Debug/`, along with a `.yak` package.

## License

[MIT](LICENSE)
