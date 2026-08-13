"""QGIS entry point for the VHRHarmonize plugin."""


def classFactory(iface):
    """Create the plugin instance expected by QGIS."""
    from .plugin import VHRHarmonizePlugin

    return VHRHarmonizePlugin(iface)
