#! /usr/bin/env python
import yaml

config_file = "config_cluster.yml"
with open(config_file, "r") as ymlfile:  # script_name without .py
    cfg = yaml.safe_load(ymlfile)
