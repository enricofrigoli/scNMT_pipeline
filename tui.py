# ATTENTION: this code is vide-coded and not thoroughly tested.
"""Textual editor for `snakeconfig.yaml`."""

from __future__ import annotations

import copy
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml
from textual.app import App, ComposeResult
from textual.containers import Horizontal, VerticalScroll
from textual.widgets import Button, Input, Select, Static, TextArea


BASE_DIR = Path(__file__).resolve().parent
CONFIG_PATH = BASE_DIR / "snakeconfig.yaml"
PIPELINES = ("star_umite", "salmon", "biscuit_methscan")


class FieldValidationError(Exception):
    """Raised when a single form field is invalid."""


@dataclass(frozen=True)
class FieldSpec:
    key: str
    label: str
    placeholder: str
    kind: str
    section: str
    must_exist: bool = False


FIELD_SPECS = (
    FieldSpec("dataset", "Dataset", "dataset name", "text", "Core"),
    FieldSpec("outdir", "Output directory", "/abs/output", "path", "Core"),
    FieldSpec(
        "pipeline",
        "Pipeline",
        "star_umite, salmon or biscuit_methscan",
        "choice",
        "Core",
    ),
    FieldSpec(
        "samples",
        "Sample filters (optional)",
        "One regex per line; leave blank to include all samples.",
        "regex_list",
        "Core",
    ),
    FieldSpec(
        "reference.genome",
        "Reference genome",
        "/abs/genome.fa",
        "path",
        "Reference",
        must_exist=True,
    ),
    FieldSpec(
        "reference.transcriptome",
        "Reference transcriptome",
        "/abs/transcriptome.fa",
        "optional_path",
        "Reference",
        must_exist=True,
    ),
    FieldSpec(
        "reference.genes",
        "Reference genes",
        "/abs/genes.gtf",
        "path",
        "Reference",
        must_exist=True,
    ),
    FieldSpec(
        "ilse_info.metadata",
        "Metadata files",
        "/abs/meta1.xls /abs/meta2.csv",
        "path_list",
        "ILSe info",
        must_exist=True,
    ),
    FieldSpec(
        "ilse_info.fastqdir",
        "FASTQ directories",
        "/abs/fastqs1 /abs/fastqs2",
        "path_list",
        "ILSe info",
        must_exist=True,
    ),
)


def load_config() -> tuple[dict[str, Any], str | None]:
    """Load an existing config file and return a warning message on failure."""
    if not CONFIG_PATH.exists():
        return {}, None

    try:
        with CONFIG_PATH.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
    except yaml.YAMLError as exc:
        return {}, f"Could not parse existing snakeconfig.yaml: {exc}"

    if not isinstance(data, dict):
        return {}, "Existing snakeconfig.yaml is not a mapping; starting blank."
    return data, None


def get_nested(
    mapping: dict[str, Any],
    dotted_key: str,
    default: str = "",
    *,
    list_separator: str = " ",
) -> str:
    """Return a nested value as a string for form prefill."""
    current: Any = mapping
    for part in dotted_key.split("."):
        if not isinstance(current, dict):
            return default
        current = current.get(part)
    if current is None:
        return default
    if isinstance(current, list):
        return list_separator.join(str(item) for item in current)
    return str(current)


def build_nested_config(base: dict[str, Any], values: dict[str, Any]) -> dict[str, Any]:
    """Merge validated UI values into a config dict."""
    config = copy.deepcopy(base)
    config["dataset"] = values["dataset"]
    config["outdir"] = values["outdir"]
    config["pipeline"] = values["pipeline"]
    config["samples"] = values["samples"]

    reference = dict(config.get("reference") or {})
    reference["genome"] = values["reference.genome"]
    reference["genes"] = values["reference.genes"]
    transcriptome = values.get("reference.transcriptome", "")
    if transcriptome:
        reference["transcriptome"] = transcriptome
    else:
        reference.pop("transcriptome", None)
    config["reference"] = reference

    ilse_info = dict(config.get("ilse_info") or {})
    ilse_info["metadata"] = " ".join(values["ilse_info.metadata"])
    ilse_info["fastqdir"] = " ".join(values["ilse_info.fastqdir"])
    config["ilse_info"] = ilse_info

    return config


def split_space_separated_paths(raw: str) -> list[str]:
    """Split a whitespace-separated path field."""
    return [part for part in raw.split() if part]


def widget_id(config_key: str) -> str:
    """Convert a config key into a Textual-safe widget id."""
    return config_key.replace(".", "-")


def validate_required_text(value: str, label: str) -> str:
    if not value.strip():
        raise FieldValidationError(f"{label} is required.")
    return value.strip()


def validate_choice(value: str, label: str, choices: tuple[str, ...]) -> str:
    cleaned = value.strip()
    if cleaned not in choices:
        raise FieldValidationError(f"{label} must be one of: {', '.join(choices)}.")
    return cleaned


def validate_sample_patterns(value: str) -> list[str]:
    """Validate one regex per nonblank line without changing pattern whitespace."""
    patterns = []
    for line_number, pattern in enumerate(value.splitlines(), start=1):
        if not pattern.strip():
            continue
        try:
            re.compile(pattern)
        except re.error as exc:
            raise FieldValidationError(
                f"Invalid sample regex on line {line_number}: {exc}"
            ) from exc
        patterns.append(pattern)
    return patterns


def validate_path(value: str, label: str, *, must_exist: bool) -> str:
    cleaned = validate_required_text(value, label)
    if not os.path.isabs(cleaned):
        raise FieldValidationError(f"{label} must be an absolute path.")
    if must_exist and not os.path.exists(cleaned):
        raise FieldValidationError(f"{label} does not exist: {cleaned}")
    return cleaned


def validate_path_list(value: str, label: str, *, must_exist: bool) -> list[str]:
    paths = split_space_separated_paths(value)
    if not paths:
        raise FieldValidationError(f"{label} cannot be empty.")
    for path in paths:
        if not os.path.isabs(path):
            raise FieldValidationError(f"{label} entries must be absolute paths.")
        if must_exist and not os.path.exists(path):
            raise FieldValidationError(f"{label} path does not exist: {path}")
    return paths


class ConfigEditor(App):
    """Textual editor for the pipeline config."""

    CSS = """
    Screen {
        align: center top;
    }

    #form {
        width: 100%;
        max-width: 120;
        padding: 1 2;
    }

    .section {
        margin-top: 1;
        text-style: bold;
    }

    .row {
        height: auto;
        margin-top: 1;
    }

    .field-label {
        width: 28;
        padding-top: 1;
    }

    Input {
        width: 1fr;
    }

    TextArea {
        width: 1fr;
        height: 6;
    }

    #status {
        margin-top: 1;
        color: $text-muted;
    }

    #status.error {
        color: $error;
    }

    #buttons {
        dock: bottom;
        width: 100%;
        height: auto;
        padding: 1 2;
        background: $surface;
        content-align: center middle;
    }
    """

    def __init__(self) -> None:
        super().__init__()
        self.base_config, self.load_warning = load_config()
        self.inputs: dict[str, Any] = {}
        self.status_widget: Static | None = None

    def compose(self) -> ComposeResult:
        """Build the form from the schema."""
        with VerticalScroll(id="form"):
            yield Static("Smart-seq3 Config Editor", classes="section")
            if self.load_warning:
                self.status_widget = Static(self.load_warning, id="status", classes="error")
            else:
                self.status_widget = Static(
                    "Edit values and confirm to write snakeconfig.yaml.",
                    id="status",
                )
            yield self.status_widget

            for section in ("Core", "Reference", "ILSe info"):
                yield Static(section, classes="section")
                for spec in (item for item in FIELD_SPECS if item.section == section):
                    value = get_nested(
                        self.base_config,
                        spec.key,
                        list_separator="\n" if spec.kind == "regex_list" else " ",
                    )
                    with Horizontal(classes="row"):
                        yield Static(spec.label, classes="field-label")
                        if spec.key == "pipeline":
                            select_value = value if value in PIPELINES else "star_umite"
                            input_widget = Select.from_values(
                                PIPELINES,
                                prompt=spec.placeholder,
                                value=select_value,
                                allow_blank=False,
                                id=widget_id(spec.key),
                            )
                        elif spec.kind == "regex_list":
                            input_widget = TextArea(
                                text=value,
                                soft_wrap=False,
                                id=widget_id(spec.key),
                            )
                            input_widget.border_title = spec.placeholder
                        else:
                            input_widget = Input(
                                value=value,
                                placeholder=spec.placeholder,
                                id=widget_id(spec.key),
                            )
                        self.inputs[spec.key] = input_widget
                        yield input_widget

        with Horizontal(id="buttons"):
            yield Button("Cancel", id="cancel")
            yield Button("Confirm", id="confirm", variant="primary")

    def on_mount(self) -> None:
        """Focus the first input when the app starts."""
        first_key = FIELD_SPECS[0].key
        self.inputs[first_key].focus()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle Confirm and Cancel buttons."""
        button_id = event.button.id
        if button_id == "cancel":
            self.exit()
            return

        if button_id != "confirm":
            return

        try:
            cleaned = self.validate_inputs()
            new_config = build_nested_config(self.base_config, cleaned)
            with CONFIG_PATH.open("w", encoding="utf-8") as handle:
                yaml.safe_dump(new_config, handle, sort_keys=False)
            self.exit()
        except FieldValidationError as exc:
            self.set_status(str(exc), error=True)
        except OSError as exc:
            self.set_status(f"Failed to write snakeconfig.yaml: {exc}", error=True)

    def set_status(self, message: str, *, error: bool = False) -> None:
        """Show a message at the bottom of the form."""
        if self.status_widget is None:
            return
        self.status_widget.update(message)
        if error:
            self.status_widget.add_class("error")
        else:
            self.status_widget.remove_class("error")

    def validate_inputs(self) -> dict[str, Any]:
        """Validate fields and return cleaned values."""
        cleaned: dict[str, Any] = {}

        dataset = validate_required_text(str(self.inputs["dataset"].value), "Dataset")
        outdir = validate_path(
            str(self.inputs["outdir"].value),
            "Output directory",
            must_exist=False,
        )
        pipeline = validate_choice(
            str(self.inputs["pipeline"].value),
            "Pipeline",
            PIPELINES,
        )
        samples = validate_sample_patterns(self.inputs["samples"].text)

        genome = validate_path(
            str(self.inputs["reference.genome"].value),
            "Reference genome",
            must_exist=True,
        )
        genes = validate_path(
            str(self.inputs["reference.genes"].value),
            "Reference genes",
            must_exist=True,
        )

        transcriptome_raw = str(self.inputs["reference.transcriptome"].value)
        transcriptome = transcriptome_raw.strip()
        if pipeline == "salmon":
            transcriptome = validate_path(
                transcriptome_raw,
                "Reference transcriptome",
                must_exist=True,
            )
        elif transcriptome:
            transcriptome = validate_path(
                transcriptome_raw,
                "Reference transcriptome",
                must_exist=True,
            )

        metadata_paths = validate_path_list(
            str(self.inputs["ilse_info.metadata"].value),
            "Metadata files",
            must_exist=True,
        )
        fastqdirs = validate_path_list(
            str(self.inputs["ilse_info.fastqdir"].value),
            "FASTQ directories",
            must_exist=True,
        )

        cleaned["dataset"] = dataset
        cleaned["outdir"] = outdir
        cleaned["pipeline"] = pipeline
        cleaned["samples"] = samples
        cleaned["reference.genome"] = genome
        cleaned["reference.genes"] = genes
        cleaned["reference.transcriptome"] = transcriptome
        cleaned["ilse_info.metadata"] = metadata_paths
        cleaned["ilse_info.fastqdir"] = fastqdirs
        return cleaned


def main() -> None:
    ConfigEditor().run()


if __name__ == "__main__":
    main()
