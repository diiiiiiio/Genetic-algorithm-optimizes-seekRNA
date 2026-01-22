import json
import os

CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.json")

def load_config():
    with open(CONFIG_PATH, "r") as f:
        return json.load(f)

CONFIG = load_config()

# Global shorthand
SEQUENCES = CONFIG["sequences"]
TARGET_SITES = CONFIG["target_sites"]
MOTIF_LENGTHS = CONFIG["motif_lengths"]
GA_PARAMS = CONFIG["ga_parameters"]
NUCS = CONFIG["nucleotides"]

