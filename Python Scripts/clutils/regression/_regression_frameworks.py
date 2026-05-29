"""Shared regression framework definitions for interaction-level analyses."""

MODEL_FRAMEWORKS = {
    "nv_vs_cartctrl": {
        "conditions": ["NV", "CAR-TCTRL"],
        "terms": ["is_treated"],
        "positive_label": "CAR-TCTRL",
        "reference_label": "NV",
    },
    "nv_vs_cartec": {
        "conditions": ["NV", "CAR-TEC"],
        "terms": ["is_treated"],
        "positive_label": "CAR-TEC",
        "reference_label": "NV",
    },
    "cartctrl_vs_cartec": {
        "conditions": ["CAR-TCTRL", "CAR-TEC"],
        "terms": ["has_colitis"],
        "positive_label": "CAR-TEC",
        "reference_label": "CAR-TCTRL",
    },
    "three_group": {
        "conditions": ["NV", "CAR-TCTRL", "CAR-TEC"],
        "terms": ["is_treated", "has_colitis"],
        "positive_label": "CAR-TEC",
        "reference_label": "NV",
    },
}
