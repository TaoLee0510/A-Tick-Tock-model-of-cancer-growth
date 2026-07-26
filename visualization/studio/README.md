# ATCG3D Studio

ATCG3D Studio is the macOS desktop control surface for the 3D simulator. It
loads and edits the complete YAML tree, validates it with the C++ executable,
starts a detached run, publishes pause/resume/stop requests, monitors atomic
run status, and embeds the existing server-side ParaView/trame viewer.

The language selector in the top bar switches the complete Studio and embedded
viewer interface between English and Simplified Chinese. English is the
default for a fresh installation. The selection is saved locally and restored
when Studio is opened again; YAML schema keys remain unchanged in both
languages.

Build:

```sh
cd visualization/studio/desktop/src-tauri
cargo test --release
cargo tauri build --bundles app
```

The bundle is written to
`target/release/bundle/macos/ATCG3D Studio.app`. It contains a standard
multi-resolution macOS icon and is ad-hoc signed during the local build;
`codesign --verify --deep --strict` must pass. It is not notarized for public
distribution.

The first release is a thin desktop app. ParaView 6.1.1 and the tested trame
Python environment remain external runtimes; the app paths can be edited in
the Run panel. Closing Studio terminates the embedded viewer process and its
heartbeat, but does not terminate the detached simulator.

When the configured run directory already contains readable status, Studio
automatically starts the embedded viewer at the latest preview frame. The
preview uses a presentation-only cell radius multiplier of 4 by default so a
whole-tumour camera view remains visible. At a non-keyframe time the backend
may then reconstruct the exact full frame from schema-8 slot-journal
checkpoints; the preview remains the initial display while that server-side
upgrade is running.

On macOS the bundle permits local HTTP content only inside its WebView so the
loopback trame server can be embedded. Studio waits for the local listener
before assigning the iframe and reports an explicit loading/failure state.
Viewer diagnostics are appended to `RUN_DIRECTORY/control/viewer.log`.

The simulator writes atomic status and request files under
`RUN_DIRECTORY/control`. While the embedded viewer is alive, Studio refreshes
`viewer.attached` once per second. The simulator then overwrites one sampled
`viz/live/current.vtkhdf` and one `viz/live/vessels.vtkhdf` at the configured
wall-clock interval. These live files never accumulate and do not alter the
archived preview/full/checkpoint sampling schedule.
