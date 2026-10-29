"""
Program developed by gqtea group.
"""

import sys

# Re-entrant plot-viewer dispatch. In a packaged (frozen) build there is no
# separate ``python.exe`` to launch ``plotViewer.py`` with, so this same
# executable is re-launched with ``--plot-viewer <manifest.json>`` to draw the
# interactive figures (see DisplayPlots.display_plots). Handle that mode here,
# before the heavy Toga/OpenGL imports, so a viewer child never loads the GUI.
if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[1] == "--plot-viewer":
    import plotViewer
    plotViewer.main(sys.argv[2])
    sys.exit(0)

import asyncio
import inspect
import traceback

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER, LEFT

# The tools are NOT imported here: importing all of them (numba, OpenGL, scipy,
# matplotlib) delayed the main window by 2.5 s warm and 11 s on a cold disk
# cache. toolLoader imports each one in the background, or when its button is
# pressed; see toolLoader.py.
from toolLoader import TOOL_TABS, ToolLoader, tool_modules
from scrollPosition import keep_scroll_position_on_relayout

class gqteaWin(toga.App):
    def startup(self):
        """Construct and show the Toga application."""
        # Scrolled tool pages must not jump to the top whenever a widget changes.
        keep_scroll_position_on_relayout()

        self.tool_loader = ToolLoader(tool_modules(), on_loaded=self._tool_loaded_from_thread)
        self._opening = []  # labels of tools whose module is still being loaded

        # 1. Main outer container
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=10))

        # 2. Header Section
        header_box = toga.Box(style=Pack(direction=COLUMN, align_items=CENTER, margin_bottom=15))
        title_label = toga.Label(
            "gQTEA-1.6.0 Molecular Analysis Toolkit",
            style=Pack(font_size=18, font_weight='bold', margin_bottom=5)
        )
        welcome_label = toga.Label(
            "Select an analysis category below to begin.",
            style=Pack(font_size=12, color="#555555")
        )
        header_box.add(title_label)
        header_box.add(welcome_label)

        # Tool loading status: a spinner and one line, refreshed as the
        # background loader finishes each tool module.
        status_row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin_top=6))
        self.loading_indicator = toga.ActivityIndicator(style=Pack(margin_right=6))
        self.status_label = toga.Label("", style=Pack(font_size=10, color="#555555"))
        status_row.add(self.loading_indicator)
        status_row.add(self.status_label)
        header_box.add(status_row)
        main_box.add(header_box)

        # 3. Create the Tabbed Container
        self.tabs = toga.OptionContainer(style=Pack(flex=1, margin_top=10))

        # Helper function to generate a scrollable list of buttons for each tab
        def create_tool_tab(tool_specs):
            tab_box = toga.Box(style=Pack(direction=COLUMN, margin=10))
            for spec in tool_specs:
                btn = toga.Button(
                    spec.label,
                    on_press=self.tool_handler(spec),
                    style=Pack(margin_bottom=1, font_size=11)
                )
                tab_box.add(btn)

            return toga.ScrollContainer(content=tab_box, style=Pack(flex=1))

        # The tabs and their buttons are defined in toolLoader.TOOL_TABS.
        for title, tool_specs in TOOL_TABS:
            self.tabs.content.append(title, create_tool_tab(tool_specs))

        # 4. Contributors Section
        contributors_box = toga.Box(style=Pack(direction=COLUMN, align_items=CENTER, margin_top=15))
        
        contributors_text = (
            "  gqteaWinToga CORE DEVELOPMENT TEAM\n"
            "  Ademir J. CAMARGO      ajc@ueg.br\n"
            "  Valter H. C. SILVA     fatioleg@ueg.br\n"
            "  Solemar S. OLIVEIRA    solemar@ueg.br\n"
            "  Hamilton B. NAPOLITANO hamilton@ueg.br\n"
            "  Luciano RIBEIRO        lribeiro@ueg.br\n"
            "  Flávio O. Sanches      flavio.neto@ifg.edu.br"
        )
        
        # Using a monospace font so the columns (names and emails) align perfectly
        self.contributors_label = toga.Label(
            contributors_text,
            style=Pack(font_size=11, font_family="monospace")
        )
        
        contributors_box.add(self.contributors_label)

        # Assemble the final view
        main_box.add(self.tabs)
        main_box.add(contributors_box)

        # Main Window setup
        self.main_window = toga.MainWindow(title="gQTEA - grupo de Química Teórica e Estrutural de Anápolis")
        self.main_window.content = main_box
        self.main_window.show()
        self.update_load_status()

    async def on_running(self):
        # Preload the tools only now, so importing them never delays the
        # window's first paint. A pressed button still jumps the queue.
        self.tool_loader.start_preload()
        self.update_load_status()

    def on_exit(self):
        self.tool_loader.stop()
        return True

    # -- opening tools -------------------------------------------------------
    def tool_handler(self, spec):
        """The on_press handler for a menu button that opens ``spec``."""
        async def open_this_tool(button, **kwargs):
            await self.open_tool(spec, button, **kwargs)
        return open_this_tool

    async def open_tool(self, spec, button, **kwargs):
        """Open a tool, importing its module first if the loader has not yet.

        The tool class is called with the Button, exactly as when the class
        itself was the on_press handler.
        """
        if not self.tool_loader.is_loaded(spec.module):
            button.enabled = False
            self._opening.append(spec.label)
            self.update_load_status()
            failure = None
            try:
                await asyncio.wrap_future(self.tool_loader.request(spec.module))
            except Exception as error:
                failure = error
            finally:
                self._opening.remove(spec.label)
                button.enabled = True
                self.update_load_status()
            if failure is not None:
                traceback.print_exception(failure)
                await self.main_window.dialog(toga.ErrorDialog(
                    "Tool could not be loaded",
                    f"{spec.label} could not be loaded:\n\n{type(failure).__name__}: {failure}",
                ))
                return
        tool = getattr(self.tool_loader.request(spec.module).result(), spec.attribute)
        result = tool(button, **kwargs)
        if inspect.isawaitable(result):
            await result

    # -- loading status ------------------------------------------------------
    def _tool_loaded_from_thread(self, name, error):
        # Called on the loader thread: touch widgets only from the Toga loop.
        try:
            self.loop.call_soon_threadsafe(self.update_load_status)
        except RuntimeError:
            pass  # the loop has closed: the app is exiting

    def load_status(self):
        """(text, busy) for the status line."""
        loaded = self.tool_loader.loaded_count()
        failed = self.tool_loader.failed_count()
        total = self.tool_loader.total_count()
        if self._opening:
            return f"Opening {self._opening[-1]}...", True
        if loaded == total:
            return f"All {total} analysis tools are ready.", False
        if loaded + failed == total:
            return (f"{loaded} of {total} analysis tools are ready; {failed} could not be "
                    "loaded (opening one shows why).", False)
        return (f"Loading analysis tools in the background: {loaded} of {total} ready. "
                "You can already open any tool.", True)

    def update_load_status(self):
        text, busy = self.load_status()
        # A Label.text write re-lays out the window: write only on change.
        if self.status_label.text != text:
            self.status_label.text = text
        if busy and not self.loading_indicator.is_running:
            self.loading_indicator.start()
        elif not busy and self.loading_indicator.is_running:
            self.loading_indicator.stop()

def main():
    return gqteaWin("gQTEA Molecular Analysis Toolkit", "br.ueg.gqtea")

if __name__ == "__main__":
    app = main()
    app.main_loop()

