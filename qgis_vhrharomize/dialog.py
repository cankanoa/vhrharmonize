"""Popup interface for editing configs and running VHRHarmonize commands."""

from __future__ import annotations

import json
import shlex
import sys
from pathlib import Path

from qgis.PyQt.QtCore import QTimer, Qt, QUrl
from qgis.PyQt.QtGui import QColor, QDesktopServices, QFontDatabase, QIcon, QPalette, QTextCursor
from qgis.PyQt.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QStyle,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)
from qgis.core import QgsApplication, QgsSettings

try:
    from qgis.PyQt.Qsci import QsciLexerYAML, QsciScintilla
except ImportError:  # Some minimal QGIS packages omit QScintilla.
    QsciLexerYAML = None
    QsciScintilla = None

from .config_manager import ConfigManager
from .interpreter import python_command
from .terminal import TerminalInput, TerminalSession


def _refresh_button(parent: QWidget, tooltip: str) -> QToolButton:
    button = QToolButton(parent)
    button.setToolTip(tooltip)
    button.setAccessibleName(tooltip)
    icon = QIcon.fromTheme("view-refresh")
    if icon.isNull():
        standard_pixmap = getattr(QStyle, "SP_BrowserReload", None)
        if standard_pixmap is None:
            standard_pixmap = QStyle.StandardPixmap.SP_BrowserReload
        icon = button.style().standardIcon(standard_pixmap)
    button.setIcon(icon)
    if icon.isNull():
        button.setText("↻")
    return button


def _play_button(parent: QWidget, tooltip: str) -> QToolButton:
    button = QToolButton(parent)
    button.setToolTip(tooltip)
    button.setAccessibleName(tooltip)
    icon = QIcon.fromTheme("media-playback-start")
    if not icon.isNull():
        button.setIcon(icon)
    else:
        button.setText("▶")
    return button


def _plain_config_editor(parent: QWidget):
    fixed_font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
    editor = QPlainTextEdit(parent)
    editor.setFont(fixed_font)
    editor.setLineWrapMode(QPlainTextEdit.NoWrap)
    return editor, None


def _yaml_editor(parent: QWidget):
    """Create a YAML-aware editor, falling back to a plain Qt editor."""
    if QsciScintilla is None or QsciLexerYAML is None:
        return _plain_config_editor(parent)

    fixed_font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
    editor = QsciScintilla(parent)
    editor.setUtf8(True)
    editor.setFont(fixed_font)
    editor.setMarginsFont(fixed_font)
    editor.setMarginLineNumbers(0, True)
    editor.setMarginWidth(0, "0000")
    editor.setIndentationsUseTabs(False)
    editor.setTabWidth(2)
    editor.setAutoIndent(True)
    editor.setIndentationGuides(True)
    editor.setBraceMatching(QsciScintilla.SloppyBraceMatch)

    palette = editor.palette()
    background = palette.color(QPalette.Base)
    foreground = palette.color(QPalette.Text)
    dark = background.lightness() < 128
    colors = {
        QsciLexerYAML.Comment: QColor("#7FBA7A" if dark else "#567A55"),
        QsciLexerYAML.Identifier: QColor("#62C8FF" if dark else "#075E91"),
        QsciLexerYAML.Keyword: QColor("#DFA5FF" if dark else "#7C3AED"),
        QsciLexerYAML.Number: QColor("#F4C97A" if dark else "#A35400"),
        QsciLexerYAML.Reference: QColor("#55D6BE" if dark else "#087F6B"),
        QsciLexerYAML.DocumentDelimiter: QColor("#FFD166" if dark else "#9A6700"),
        QsciLexerYAML.TextBlockMarker: QColor("#FFD166" if dark else "#9A6700"),
        QsciLexerYAML.Operator: foreground,
        QsciLexerYAML.SyntaxErrorMarker: QColor("#FFFFFF"),
    }
    styles = (
        QsciLexerYAML.Default,
        QsciLexerYAML.Comment,
        QsciLexerYAML.Identifier,
        QsciLexerYAML.Keyword,
        QsciLexerYAML.Number,
        QsciLexerYAML.Reference,
        QsciLexerYAML.DocumentDelimiter,
        QsciLexerYAML.TextBlockMarker,
        QsciLexerYAML.Operator,
        QsciLexerYAML.SyntaxErrorMarker,
    )
    lexer = QsciLexerYAML(editor)
    lexer.setDefaultColor(foreground)
    lexer.setDefaultPaper(background)
    lexer.setDefaultFont(fixed_font)
    for style in styles:
        lexer.setPaper(background, style)
        lexer.setFont(fixed_font, style)
    for style, color in colors.items():
        lexer.setColor(color, style)
    lexer.setPaper(QColor("#C42B1C"), QsciLexerYAML.SyntaxErrorMarker)
    editor.setLexer(lexer)
    editor.setCaretForegroundColor(foreground)
    editor.setMarginsForegroundColor(palette.color(QPalette.Disabled, QPalette.Text))
    editor.setMarginsBackgroundColor(background)
    return editor, lexer


class ConfigEditor(QWidget):
    """A text editor which saves to its config file after a short debounce."""

    def __init__(
        self,
        manager: ConfigManager,
        filename: str,
        *,
        reset_at_right: bool,
        parent: QWidget | None = None,
    ):
        super().__init__(parent)
        self.manager = manager
        self.filename = filename
        editor_factory = (
            _yaml_editor if filename.lower().endswith((".yml", ".yaml")) else _plain_config_editor
        )
        self.editor, self._lexer = editor_factory(self)
        self._set_text(manager.read(filename))
        self.editor.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self._save_timer = QTimer(self)
        self._save_timer.setSingleShot(True)
        self._save_timer.setInterval(400)
        self._save_timer.timeout.connect(self.flush)
        self.editor.textChanged.connect(self._save_timer.start)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.editor, 1)
        if reset_at_right:
            reset = _refresh_button(self, f"Reset {filename} to its default")
            reset.clicked.connect(self.reset)
            layout.addWidget(reset, 0, Qt.AlignTop)

    def flush(self) -> None:
        self._save_timer.stop()
        text = self._text()
        if text != self.manager.read(self.filename):
            self.manager.write(self.filename, text)

    def reset(self) -> None:
        self._save_timer.stop()
        self._set_text(self.manager.reset(self.filename))

    def _text(self) -> str:
        if self._lexer is not None:
            return self.editor.text()
        return self.editor.toPlainText()

    def _set_text(self, text: str) -> None:
        if self._lexer is not None:
            self.editor.setText(text)
        else:
            self.editor.setPlainText(text)


class VHRHarmonizeDialog(QDialog):
    """Main VHRHarmonize interface."""

    SETTINGS_PROVIDER = "vhrharmonize/provider"

    def __init__(self, plugin_dir: Path, manager: ConfigManager, parent: QWidget | None = None):
        super().__init__(parent)
        self.plugin_dir = Path(plugin_dir)
        self.manager = manager
        self.settings = QgsSettings()
        self.interpreter = python_command()
        py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        self.dependency_target = (
            Path(QgsApplication.qgisSettingsDirPath()) / "python" / "dependencies" / py_version
        )
        self.terminal = TerminalSession(self.plugin_dir, self)
        self.optional_command_inputs: dict[str, QLineEdit] = {}
        self.optional_command_defaults: dict[str, str] = {}

        self.setWindowTitle("VHRHarmonize")
        self.resize(980, 720)
        root = QVBoxLayout(self)
        self.tabs = QTabWidget(self)
        root.addWidget(self.tabs)

        self.setup_tab = self._build_setup_tab()
        self.provider_tab = self._build_provider_tab()
        self.hpc_tab = self._build_hpc_tab()
        self.run_tab = self._build_run_tab()
        self.tabs.addTab(self.setup_tab, "Setup")
        self.tabs.addTab(self.provider_tab, "Provider")
        self.tabs.addTab(self.hpc_tab, "High Performance Compute")
        self.tabs.addTab(self.run_tab, "Run")

        self.terminal.output_received.connect(self._append_output)
        self.terminal.input_echo_changed.connect(self._set_terminal_echo)
        self.terminal.state_changed.connect(self.terminal_status.setText)
        self.terminal.start()

    def _build_setup_tab(self) -> QWidget:
        page = QWidget(self)
        layout = QVBoxLayout(page)
        layout.addWidget(QLabel(f"QGIS Python interpreter: {self.interpreter}", page))
        layout.addWidget(QLabel(f"qpip dependency directory: {self.dependency_target}", page))
        explanation = QLabel(
            "Optional dependency commands are editable. Press play to place a command in the "
            "interactive terminal input on the Run tab; press Enter there to run it.",
            page,
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)

        scroll = QScrollArea(page)
        scroll.setWidgetResizable(True)
        content = QWidget(scroll)
        content_layout = QVBoxLayout(content)
        manifest_path = self.plugin_dir / "optional_dependency_commands.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for group in manifest.get("groups", []):
            name = str(group["name"])
            tokens = [
                str(token).format(
                    python=self.interpreter,
                    target=str(self.dependency_target),
                )
                for token in group["command"]
            ]
            default_command = shlex.join(tokens)
            self.optional_command_defaults[name] = default_command

            section = QWidget(content)
            section_layout = QVBoxLayout(section)
            section_layout.setContentsMargins(0, 4, 0, 4)
            header = QHBoxLayout()
            name_label = QLabel(name, section)
            name_label.setStyleSheet("font-weight: bold")
            header.addWidget(name_label)
            reset = _refresh_button(section, f"Reset the {name} install command")
            reset.clicked.connect(
                lambda _checked=False, group_name=name: self._reset_optional_command(group_name)
            )
            header.addWidget(reset)
            header.addStretch(1)
            section_layout.addLayout(header)

            command_row = QHBoxLayout()
            command_input = QLineEdit(default_command, section)
            command_input.setToolTip(f"Editable pip command for the {name} optional dependency group")
            self.optional_command_inputs[name] = command_input
            command_row.addWidget(command_input, 1)
            play = _play_button(section, f"Open the {name} install command in the terminal")
            play.clicked.connect(
                lambda _checked=False, group_name=name: self._load_optional_dependency(group_name)
            )
            command_row.addWidget(play)
            section_layout.addLayout(command_row)
            content_layout.addWidget(section)

        content_layout.addStretch(1)
        scroll.setWidget(content)
        layout.addWidget(scroll, 1)
        self.setup_status = QLabel("", page)
        layout.addWidget(self.setup_status)
        return page

    def _build_provider_tab(self) -> QWidget:
        page = QWidget(self)
        layout = QVBoxLayout(page)
        options = QHBoxLayout()
        options.addWidget(QLabel("Provider:", page))
        self.provider_combo = QComboBox(page)
        self.provider_combo.addItem("WorldView", "worldview")
        saved = str(self.settings.value(self.SETTINGS_PROVIDER, "WorldView"))
        index = self.provider_combo.findText(saved)
        self.provider_combo.setCurrentIndex(index if index >= 0 else 0)
        self.settings.setValue(self.SETTINGS_PROVIDER, self.provider_combo.currentText())
        self.provider_combo.currentTextChanged.connect(self._provider_changed)
        options.addWidget(self.provider_combo)
        reset = _refresh_button(page, "Reset the provider config to its default")
        reset.clicked.connect(lambda: self.provider_editor.reset())
        options.addWidget(reset)
        options.addStretch(1)
        layout.addLayout(options)

        self.provider_editor = ConfigEditor(
            self.manager,
            ConfigManager.WORLDVIEW,
            reset_at_right=False,
            parent=page,
        )
        layout.addWidget(self.provider_editor, 1)
        return page

    def _build_hpc_tab(self) -> QWidget:
        page = QWidget(self)
        layout = QVBoxLayout(page)
        inner_tabs = QTabWidget(page)
        self.hpc_editor = ConfigEditor(
            self.manager,
            ConfigManager.HPC,
            reset_at_right=True,
            parent=inner_tabs,
        )
        self.slurm_editor = ConfigEditor(
            self.manager,
            ConfigManager.SLURM,
            reset_at_right=True,
            parent=inner_tabs,
        )
        inner_tabs.addTab(self.hpc_editor, "HPC")
        inner_tabs.addTab(self.slurm_editor, "Slurm")
        layout.addWidget(inner_tabs)
        return page

    def _build_run_tab(self) -> QWidget:
        page = QWidget(self)
        layout = QVBoxLayout(page)
        self.run_hpc = QCheckBox("Run Pipeline High Performance Compute", page)
        self.run_hpc.setChecked(False)
        layout.addWidget(self.run_hpc)

        self.run_modes = QStackedWidget(page)
        self.run_modes.addWidget(self._build_local_controls())
        self.run_modes.addWidget(self._build_hpc_controls())
        self.run_hpc.toggled.connect(lambda checked: self.run_modes.setCurrentIndex(1 if checked else 0))
        layout.addWidget(self.run_modes)

        self.output_log = QPlainTextEdit(page)
        self.output_log.setReadOnly(True)
        self.output_log.setPlaceholderText("The interactive shell will appear here.")
        self.output_log.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))

        status_row = QHBoxLayout()
        self.terminal_status = QLabel("Starting interactive shell…", page)
        status_row.addWidget(self.terminal_status)
        status_row.addStretch(1)
        clear_output = QPushButton("Clear Output", page)
        clear_output.clicked.connect(self.output_log.clear)
        status_row.addWidget(clear_output)
        restart = QPushButton("Restart Shell", page)
        restart.setToolTip("Start a fresh shell session")
        restart.clicked.connect(self.terminal.restart)
        status_row.addWidget(restart)
        layout.addLayout(status_row)

        layout.addWidget(self.output_log, 1)

        terminal_row = QHBoxLayout()
        terminal_row.addWidget(QLabel("$", page))
        self.terminal_input = TerminalInput(page)
        self.terminal_input.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        self.terminal_input.setPlaceholderText(
            "Type a command or response; Enter sends it directly to the shell"
        )
        self.terminal_input.submitted.connect(self._send_terminal_line)
        terminal_row.addWidget(self.terminal_input, 1)
        send = QPushButton("Send", page)
        send.setToolTip("Send the input directly to the interactive shell")
        send.clicked.connect(self.terminal_input.submit)
        terminal_row.addWidget(send)
        interrupt = QPushButton("Ctrl+C", page)
        interrupt.setToolTip("Interrupt the active terminal command")
        interrupt.clicked.connect(self.terminal.interrupt)
        terminal_row.addWidget(interrupt)
        eof = QPushButton("Ctrl+D", page)
        eof.setToolTip("Send end-of-file to the interactive shell")
        eof.clicked.connect(self.terminal.send_eof)
        terminal_row.addWidget(eof)
        layout.addLayout(terminal_row)
        return page

    def _build_local_controls(self) -> QWidget:
        page = QWidget(self)
        layout = QGridLayout(page)
        start = QPushButton("Start", page)
        layout.setRowMinimumHeight(0, start.sizeHint().height())
        start.setToolTip("Load the provider pipeline command into the command field")
        start.clicked.connect(self._set_local_command_preset)
        layout.addWidget(start, 1, 0)
        open_configs = QPushButton("Open Configs Folder", page)
        open_configs.setToolTip("Open the editable config directory in the system file manager")
        open_configs.clicked.connect(self._open_configs_folder)
        layout.addWidget(open_configs, 1, 1)
        layout.setColumnStretch(2, 1)
        layout.setColumnStretch(3, 1)
        return page

    def _build_hpc_controls(self) -> QWidget:
        page = QWidget(self)
        layout = QGridLayout(page)
        layout.addWidget(QLabel("Staged HPC file:", page), 0, 0)
        self.staged_hpc_input = QLineEdit(
            self._display_staged_hpc_path(self.manager.staged_hpc_path()),
            page,
        )
        self.staged_hpc_input.setPlaceholderText("1.staged.hpc.yml")
        self.staged_hpc_input.setToolTip(
            "Filename used by HPC commands after Prepare. Relative names are resolved in configs/."
        )
        layout.addWidget(self.staged_hpc_input, 0, 1, 1, 2)
        open_configs = QPushButton("Open Configs Folder", page)
        open_configs.setToolTip("Open the editable config directory in the system file manager")
        open_configs.clicked.connect(self._open_configs_folder)
        layout.addWidget(open_configs, 0, 3)

        actions = (
            ("Prepare", "prepare"),
            ("Upload", "upload"),
            ("Start", "start"),
            ("Status", "status"),
            ("Download", "download"),
            ("Stop Job", "stop"),
            ("Close Connection", "close"),
        )
        for index, (label, action) in enumerate(actions):
            button = QPushButton(label, page)
            button.clicked.connect(lambda _checked=False, value=action: self._set_hpc_command_preset(value))
            layout.addWidget(button, 1 + index // 4, index % 4)
        full_run = QPushButton("Full Run", page)
        full_run.setToolTip("Load prepare, upload, start, and status into the command field")
        full_run.clicked.connect(self._set_full_run_preset)
        layout.addWidget(full_run, 2, 3)
        return page

    def _open_configs_folder(self) -> None:
        folder_url = QUrl.fromLocalFile(str(self.manager.config_dir.resolve()))
        if not QDesktopServices.openUrl(folder_url):
            self.terminal_status.setText(
                f"Could not open the configs folder: {self.manager.config_dir}"
            )

    def _provider_changed(self, provider: str) -> None:
        self.settings.setValue(self.SETTINGS_PROVIDER, provider)

    def _reset_optional_command(self, name: str) -> None:
        self.optional_command_inputs[name].setText(self.optional_command_defaults[name])

    def _load_optional_dependency(self, name: str) -> None:
        text = self.optional_command_inputs[name].text().strip()
        if not text:
            self.setup_status.setText(f"The {name} command is empty.")
            return
        self._set_terminal_text(text, switch_to_run=True)
        self.setup_status.setText(f"Loaded {name} in the terminal. Press Enter to install it.")

    def _flush_configs(self) -> None:
        self.provider_editor.flush()
        self.hpc_editor.flush()
        self.slurm_editor.flush()

    def _module_command(self, module: str, args: list[str]) -> str:
        return shlex.join([self.interpreter, "-u", "-m", module, *args])

    def _local_command(self) -> str:
        config = str(self.manager.path(ConfigManager.WORLDVIEW))
        return self._module_command(
            "vhrharmonize.cli.worldview",
            ["--config-yaml", config],
        )

    def _hpc_command(self, action: str) -> str:
        config_path = (
            self.manager.path(ConfigManager.HPC)
            if action == "prepare"
            else self._selected_staged_hpc_path()
        )
        return self._module_command(
            "vhrharmonize.slurm",
            [action, "--config", str(config_path)],
        )

    def _set_local_command_preset(self) -> None:
        self._set_terminal_text(self._local_command())

    def _set_hpc_command_preset(self, action: str) -> None:
        if action == "prepare":
            self._refresh_staged_hpc_file()
        self._set_terminal_text(self._hpc_command(action))

    def _set_full_run_preset(self) -> None:
        self._refresh_staged_hpc_file()
        commands = [
            self._hpc_command(action)
            for action in ("prepare", "upload", "start", "status")
        ]
        self._set_terminal_text(" && ".join(commands))

    def _refresh_staged_hpc_file(self) -> Path:
        """Update the staged filename from the output declared by the HPC config."""
        self.hpc_editor.flush()
        path = self.manager.staged_hpc_path()
        self.staged_hpc_input.setText(self._display_staged_hpc_path(path))
        return path

    def _display_staged_hpc_path(self, path: Path) -> str:
        try:
            return path.relative_to(self.manager.config_dir).as_posix()
        except ValueError:
            pass
        try:
            return path.relative_to(self.plugin_dir).as_posix()
        except ValueError:
            return str(path)

    def _selected_staged_hpc_path(self) -> Path:
        text = self.staged_hpc_input.text().strip()
        if not text:
            return self.manager.staged_hpc_path()
        path = Path(text).expanduser()
        if path.is_absolute():
            return path
        if path.parts and path.parts[0] == self.manager.config_dir.name:
            return self.plugin_dir / path
        return self.manager.config_dir / path

    def _set_terminal_text(self, text: str, *, switch_to_run: bool = False) -> None:
        self.terminal_input.setText(text)
        self.terminal_input.setCursorPosition(len(text))
        if switch_to_run:
            self.tabs.setCurrentWidget(self.run_tab)
        self.terminal_input.setFocus()

    def _send_terminal_line(self, text: str) -> None:
        self._flush_configs()
        self.terminal.send_line(text)

    def _set_terminal_echo(self, enabled: bool) -> None:
        normal_mode = getattr(QLineEdit, "Normal", None)
        password_mode = getattr(QLineEdit, "Password", None)
        if normal_mode is None:
            normal_mode = QLineEdit.EchoMode.Normal
            password_mode = QLineEdit.EchoMode.Password
        self.terminal_input.setEchoMode(normal_mode if enabled else password_mode)
        self.terminal_input.setPlaceholderText(
            "Type a command or response; Enter sends it directly to the shell"
            if enabled
            else "Password/input — terminal echo is disabled"
        )

    def _append_output(self, text: str) -> None:
        cursor = self.output_log.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text)
        self.output_log.setTextCursor(cursor)
        self.output_log.ensureCursorVisible()

    def closeEvent(self, event) -> None:
        self._flush_configs()
        super().closeEvent(event)
