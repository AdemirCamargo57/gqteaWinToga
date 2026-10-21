"""Keep a Toga ScrollContainer's scroll position across re-layouts (WinForms).

toga-winforms 0.5's ``ScrollContainer.refreshed`` ends with
``AutoScroll = False; AutoScroll = True``, a workaround for a horizontal scroll
bar that stayed visible after a DPI change. In WinForms turning AutoScroll off
scrolls the panel back to the top at once (its pixels are moved on screen
immediately), and ``refreshed`` runs on every re-layout of the page -- which
Toga does whenever a widget inside changes size-relevant state, such as a
Label's text or a Selection's items. So choosing an option far down a scrolled
tab (the 3D viewer's QTAIM Canvas section, for one) made the page flash and
land at the top. Restoring the position afterwards does not help: right after
the toggle the panel has no scroll range yet, so WinForms clamps it to 0.

``keep_scroll_position_on_relayout`` replaces ``refreshed`` with the same
layout pass minus the toggle. Toga's original, toggle included, still runs
whenever the DPI scale has changed since the container's last layout (the case
the toggle exists for), followed by a restore of the position.

The replacement mirrors toga-winforms 0.5's code, so it is installed only for
that release line. Toga is imported only inside the function, so importing
this module stays GUI-free.
"""

SUPPORTED_TOGA_WINFORMS = "0.5."


def keep_scroll_position_on_relayout() -> bool:
    """Patch toga-winforms' ScrollContainer once; True when the patch is in place.

    Returns False (and changes nothing) where the WinForms backend is missing
    or is not a 0.5 release. Safe to call more than once.
    """
    try:
        import toga_winforms
        from System.Drawing import Point
        from System.Windows.Forms import SystemInformation
        from toga_winforms.widgets.scrollcontainer import ScrollContainer
        from travertino.node import Node
    except ImportError:
        return False
    if not str(getattr(toga_winforms, "__version__", "")).startswith(SUPPORTED_TOGA_WINFORMS):
        return False

    original = ScrollContainer.refreshed
    if getattr(original, "keeps_scroll_position", False):
        return True

    def refreshed(self):
        dpi_scale = self.dpi_scale
        if getattr(self, "_layout_dpi_scale", None) != dpi_scale:
            # First layout or a DPI change: Toga's own pass, toggle and all.
            self._layout_dpi_scale = dpi_scale
            position = self.native.AutoScrollPosition
            original(self)
            if position.X or position.Y:
                self.native.PerformLayout()   # rebuild the scroll range the toggle cleared
                self.native.AutoScrollPosition = Point(-position.X, -position.Y)
            return

        # toga-winforms 0.5 ScrollContainer.refreshed, without the AutoScroll toggle.
        full_width, full_height = (self.native_width, self.native_height)
        inset_width = full_width - SystemInformation.VerticalScrollBarWidth
        inset_height = full_height - SystemInformation.HorizontalScrollBarHeight
        layout = self.interface.content.layout

        def apply_insets():
            need_scrollbar = False
            if self.vertical and (layout.height > self.height):
                need_scrollbar = True
                self.native_width = inset_width
            if self.horizontal and (layout.width > self.width):
                need_scrollbar = True
                self.native_height = inset_height
            return need_scrollbar

        if apply_insets():
            Node.refresh(self.interface.content, self)
            apply_insets()

        self.apply_layout(
            layout.width if self.horizontal else 0,
            layout.height if self.vertical else 0,
        )
        self.native_width, self.native_height = full_width, full_height

    refreshed.keeps_scroll_position = True
    refreshed.__doc__ = original.__doc__
    ScrollContainer.refreshed = refreshed
    return True
