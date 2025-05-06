from shiny import App, ui, render, reactive
import pandas as pd
from collections import Counter
from scipy.stats import chisquare

# Define Mendelian models
MENDELIAN_MODELS = {
    "3:1": [3, 1],
    "1:2:1": [1, 2, 1],
    "1:1:1:1": [1, 1, 1, 1],
    "9:3:3:1": [9, 3, 3, 1],
    "12:3:1": [12, 3, 1],
    "9:7": [9, 7],
    "15:1": [15, 1],
    "9:3:4": [9, 3, 4],
    "13:3": [13, 3],
    "9:6:1": [9, 6, 1]
}

# Define UI
app_ui = ui.page_fluid(
    ui.h2("Mendelian Ratio Chi-square Tester"),
    ui.markdown("**Paste your phenotypic data (one value per line):**"),
    ui.input_text_area("phenodata", "", rows=15),
    ui.hr(),
    ui.output_ui("result_ui")
)

# Define server logic
def server(input, output, session):

    @reactive.calc
    def parsed_data():
        lines = input.phenodata().strip().splitlines()
        cleaned = [line.strip().lower() for line in lines if line.strip()]
        return Counter(cleaned)

    @reactive.calc
    def matched_models():
        count = len(parsed_data())
        return {
            name: ratio
            for name, ratio in MENDELIAN_MODELS.items()
            if len(ratio) == count
        }

    @reactive.calc
    def test_results():
        obs
