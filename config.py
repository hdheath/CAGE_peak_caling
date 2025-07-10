# config.py
import tomli  # or: import tomli as toml if you renamed it

def load_config(path):
    with open(path, "rb") as f:
        return tomli.load(f) or {}
