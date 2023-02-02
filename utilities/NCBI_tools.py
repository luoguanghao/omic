import numpy as np
import pandas as pd
from pathlib import Path
import os, sys
import re

from bs4 import BeautifulSoup
import urllib.request


def GSM2name(url):
    page = urllib.request.urlopen(url)
    soup = BeautifulSoup(page,"html.parser")
    name = str(soup).split('\n')[1].strip().split(' = ')[1]
    return name



