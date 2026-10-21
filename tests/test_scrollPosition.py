"""scrollPosition: Toga ScrollContainers keep their scroll position across re-layouts."""
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

import scrollPosition

SRC = Path(__file__).resolve().parents[1]


def test_import_does_not_load_toga():
    code = "import sys, scrollPosition; assert 'toga' not in sys.modules and 'toga_winforms' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


@pytest.fixture
def scroll_container(monkeypatch):
    """toga-winforms' ScrollContainer class, with refreshed restored after the test."""
    sc = pytest.importorskip("toga_winforms.widgets.scrollcontainer")
    monkeypatch.setattr(sc.ScrollContainer, "refreshed", sc.ScrollContainer.refreshed)
    return sc.ScrollContainer


class FakeNative:
    """The scrolling Panel: records AutoScroll writes and the position assigned."""

    def __init__(self, x=0, y=0):
        self.AutoScrollPosition = SimpleNamespace(X=-x, Y=-y)
        self.auto_scroll_writes, self.layouts = [], 0

    def __setattr__(self, name, value):
        if name == "AutoScroll":
            self.auto_scroll_writes.append(value)
        object.__setattr__(self, name, value)

    def PerformLayout(self):
        self.layouts += 1


def fake_container(native, dpi=1.0, content_height=100, height=400):
    layouts = []
    c = SimpleNamespace(
        native=native, dpi_scale=dpi, native_width=500, native_height=height, width=500, height=height,
        vertical=True, horizontal=False,
        interface=SimpleNamespace(content=SimpleNamespace(layout=SimpleNamespace(width=480, height=content_height))),
        apply_layout=lambda w, h: layouts.append((w, h)),
    )
    c.layouts = layouts
    return c


class TestPatch:
    def test_installs_once(self, scroll_container):
        assert scrollPosition.keep_scroll_position_on_relayout()
        patched = scroll_container.refreshed
        assert patched.keeps_scroll_position
        assert scrollPosition.keep_scroll_position_on_relayout() and scroll_container.refreshed is patched

    def test_other_toga_versions_are_left_alone(self, scroll_container, monkeypatch):
        import toga_winforms
        monkeypatch.setattr(toga_winforms, "__version__", "0.6.0")
        original = scroll_container.refreshed
        assert scrollPosition.keep_scroll_position_on_relayout() is False
        assert scroll_container.refreshed is original

    def test_relayout_never_toggles_auto_scroll(self, scroll_container):
        scrollPosition.keep_scroll_position_on_relayout()
        native = FakeNative(y=746)
        c = fake_container(native)
        c._layout_dpi_scale = 1.0            # laid out before at this DPI
        scroll_container.refreshed(c)
        assert native.auto_scroll_writes == []                  # the toggle is what reset the position
        assert native.AutoScrollPosition.Y == -746
        assert c.layouts == [(0, 100)]                          # the same layout pass as Toga's
        assert (c.native_width, c.native_height) == (500, 400)

    def test_tall_content_keeps_the_scroll_bar_inset(self, scroll_container, monkeypatch):
        scrollPosition.keep_scroll_position_on_relayout()
        import travertino.node
        calls = []
        monkeypatch.setattr(travertino.node.Node, "refresh", lambda content, container: calls.append(
            container.native_width))
        c = fake_container(FakeNative(y=10), content_height=2000)
        c._layout_dpi_scale = 1.0
        scroll_container.refreshed(c)
        assert calls and calls[0] < 500                         # laid out narrower, beside the scroll bar
        assert c.native_width == 500                            # and the full width restored afterwards

    def test_dpi_change_runs_togas_pass_and_restores_the_position(self, scroll_container, monkeypatch):
        seen = []
        monkeypatch.setattr(scroll_container, "refreshed", lambda self: seen.append(self))
        scrollPosition.keep_scroll_position_on_relayout()
        native = FakeNative(y=300)
        c = fake_container(native, dpi=1.25)
        c._layout_dpi_scale = 1.0
        scroll_container.refreshed(c)
        assert seen == [c] and c._layout_dpi_scale == 1.25
        assert native.layouts == 1
        assert (native.AutoScrollPosition.X, native.AutoScrollPosition.Y) == (0, 300)   # Point(x, y), positive
