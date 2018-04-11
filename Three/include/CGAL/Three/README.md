# Review of Three

## Main goals of the new API
- ease the maintenance of the code,
- ease the maintenance of binary compatibility,
- modularity,
- allow to setup an API for script languages (using Qt5 MOC system).

## `Scene`, `Viewer`, `Scene_item`, `Scene_interface`

I wonder why I had decided to use that `Scene_interface` interface. The demo in itself is so heavy that decoupling the implementation of `Scene` and `Viewer` now seems to be too much.

- We want to decouple the `Scene_item` and `Viewer`, as much as possible.
- I do not think plugins should see the viewer implementation details as well. However, some plugins need to draw trans

## Plugins

The plugin system is the tool that allows the modularity of the application.

## I/O Plugins

The I/O system should be rethinked completely, now that we have real I/O formats, and real use cases. The main issue of the old API was that the `load()` function was returning a `Scene_item*`: it assumed that one file is exacly one item.
