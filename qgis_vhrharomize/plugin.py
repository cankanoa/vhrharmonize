"""QGIS plugin lifecycle implementation."""

from pathlib import Path

from qgis.core import QgsSettings
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction

from .config_manager import ConfigManager


class VHRHarmonizePlugin:
    """Register and display the VHRHarmonize popup dialog."""

    def __init__(self, iface):
        self.iface = iface
        self.plugin_dir = Path(__file__).resolve().parent
        self.config_manager = ConfigManager(self.plugin_dir)
        self.config_manager.ensure_user_configs()
        settings = QgsSettings()
        if not settings.contains("vhrharmonize/provider"):
            settings.setValue("vhrharmonize/provider", "WorldView")
        self.action = None
        self.dialog = None

    def initGui(self) -> None:
        icon = QIcon(str(self.plugin_dir / "icon.png"))
        self.action = QAction(icon, "VHRHarmonize", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("VHRHarmonize", self.action)

    def unload(self) -> None:
        if self.dialog is not None:
            self.dialog.terminal.shutdown()
            self.dialog.close()
        if self.action is not None:
            self.iface.removeToolBarIcon(self.action)
            self.iface.removePluginMenu("VHRHarmonize", self.action)

    def run(self) -> None:
        if self.dialog is None:
            from .dialog import VHRHarmonizeDialog

            self.dialog = VHRHarmonizeDialog(
                self.plugin_dir,
                self.config_manager,
                self.iface.mainWindow(),
            )
        self.dialog.show()
        self.dialog.raise_()
        self.dialog.activateWindow()
