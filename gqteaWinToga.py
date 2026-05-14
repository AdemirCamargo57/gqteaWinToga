"""
Program developed by gqtea group.
"""

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER, LEFT

# --- Existing imports ---
from bond import BondUI
from allBondAnalysis import allBondAnalysisUI
from allAnglesAnalysis import allAnglesAnalysisUI
from allDihedralAnalysis import allDihedralAnalysisUI
from hbond import HBondUI
from bondAngle import BondAngleUI
from dihedralAngle import DihedralUI
from cpmdInput import CpmdInputUI
from spectraSH import SurfaceHoppingUI
from collision import CollisionUI
from selectFrames import SelectFramesUI
from plotter import PlotterUI
from cpmdInputToXYZ import CpmdInputToXYZUI
from radialDistribution import RadialFunctionUI
from meanResidenceTime import MeanResidenceTimeUI
from legacy_meanResidenceTime import LegacyMeanResidenceTimeUI
from molecularViewer import MolecularViewerUI
from autocorrelationFunction import AutoCorrelationFunctionUI
from mixtureSolventBox import MixtureSolventBoxUI
from single_solute_solvent_box import SingleSoluteSolventBoxUI
from sharedWallBox import SharedWallBoxUI
from coordinateConverter import CoordinateConverterUI
from classicalRate import ClassicalRateUI
from cpTrajec2xyz import CPtraj2xyzUI
from cpxForces import CPforcesUI
from xyz2bohrConverter import XYZ2BohrConverterUI
from rangeFramesSelection import RangeFramesSelectionUI
from sh_geom_analyzer import SHGeomAnalyzerUI
from runatom_input_gen import RunatomGenerator
from unitConverter import UnitConvertUI
from molecularAxisAlignment import MolecularAxisAlignmentUI
from cpx_input_builder import CPInputBuilderUI
from orca_input_builder import ORCAInputBuilderUI
from gqteaMDinputBuilder import GqteaMDInputBuilderUI

class gqteaWin(toga.App):
    def startup(self):
        """Construct and show the Toga application."""
        
        # 1. Main outer container
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=10))

        # 2. Header Section
        header_box = toga.Box(style=Pack(direction=COLUMN, align_items=CENTER, margin_bottom=15))
        title_label = toga.Label(
            "gQTEA-0.3.7 Molecular Analysis Toolkit",
            style=Pack(font_size=18, font_weight='bold', margin_bottom=5)
        )
        welcome_label = toga.Label(
            "Select an analysis category below to begin.",
            style=Pack(font_size=12, color="#555555")
        )
        header_box.add(title_label)
        header_box.add(welcome_label)
        main_box.add(header_box)

        # 3. Create the Tabbed Container
        self.tabs = toga.OptionContainer(style=Pack(flex=1, margin_top=10))

        # Helper function to generate a scrollable list of buttons for each tab
        def create_tool_tab(button_definitions):
            tab_box = toga.Box(style=Pack(direction=COLUMN, margin=10))
            for label, action in button_definitions:
                btn = toga.Button(
                    label, 
                    on_press=action, 
                    style=Pack(margin_bottom=1, font_size=11)
                )
                tab_box.add(btn)
            
            return toga.ScrollContainer(content=tab_box, style=Pack(flex=1))

        # --- Define the layout for each tab ---        
        geometry_tools = [
            ("Bond length analysis", BondUI),
            ("All bond distance analysis", allBondAnalysisUI),
            ("Bond angle analysis", BondAngleUI),
            ("All bond angle analysis", allAnglesAnalysisUI),
            ("Dihedral angle analysis", DihedralUI),
            ("All dihedral angle analysis", allDihedralAnalysisUI),
            ("Hydrogen bond analysis", HBondUI),
        ]

        cpmd_tools = [
            ("CPMD Inputs", CpmdInputUI),
            ("SH Input Builder", SurfaceHoppingUI),
            ("Collision Input", CollisionUI),
            ("cp.x Input Builder", CPInputBuilderUI),
            ("ORCA Input Builder", ORCAInputBuilderUI),
            ("gqteaMD Input Builder", GqteaMDInputBuilderUI),
            ("Vanderbilt runatom.x Input Builder", RunatomGenerator),
        ]

        structural_tools = [
            ("Radial Distribution Function", RadialFunctionUI),
            ("Mean Residence Time", MeanResidenceTimeUI),
            ("Legacy Mean Residence Time", LegacyMeanResidenceTimeUI),
            ("Autocorrelation function", AutoCorrelationFunctionUI),
            ("Single Solute solvent box", SingleSoluteSolventBoxUI),
            ("Mixture of two solvent box", MixtureSolventBoxUI),
            ("Shared-wall double solvent box", SharedWallBoxUI),
        ]

        Thermo_tools = [
            ("Classical rate constant", ClassicalRateUI),
        ]

        general_tools = [
            ("3D Molecular Viewer", MolecularViewerUI),
            ("Energy plots", PlotterUI),
            ("Molecular Axis Alignment", MolecularAxisAlignmentUI),
            ("Select Frames", SelectFramesUI),
            ("Frame selection by interatomic distance range", RangeFramesSelectionUI),
            ("CPMD Input to XYZ Converter", CpmdInputToXYZUI),
            ("SH Geometry Analyzer", SHGeomAnalyzerUI),
            ("Convert *.pos file to trajec.xyz", CPtraj2xyzUI),
            ("Compute forces from cp.x *.for file", CPforcesUI),
            ("Coordinate converter", CoordinateConverterUI),
            ("Convert coordinates from Å to Bohr", XYZ2BohrConverterUI),
            ("Unit Converter", UnitConvertUI),            
        ]

        # Add the generated tabs to the OptionContainer 
        self.tabs.content.append("Geometry", create_tool_tab(geometry_tools))
        self.tabs.content.append("Inputs", create_tool_tab(cpmd_tools))
        self.tabs.content.append("Structural", create_tool_tab(structural_tools))
        self.tabs.content.append("Thermo", create_tool_tab(Thermo_tools))
        self.tabs.content.append("Tools", create_tool_tab(general_tools))

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

def main():
    return gqteaWin("gQTEA Molecular Analysis Toolkit", "br.ueg.gqtea")

if __name__ == "__main__":
    app = main()
    app.main_loop()

