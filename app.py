from shiny import App, reactive, render, ui
import pandas as pd
import numpy as np
from scipy.stats import chisquare
from plotnine import ggplot, aes, geom_bar, labs, theme_minimal, scale_fill_manual, position_dodge
from pathlib import Path

# Define segregation models and their explanations
models = {
    "3:1": ([3, 1], "Monogenic inheritance: dominant vs recessive (3 = dominant, 1 = recessive)."),
    "1:2:1": ([1, 2, 1], "Monogenic inheritance with incomplete dominance or codominance."),
    "1:1:1:1": ([1, 1, 1, 1], "Dihybrid cross with independent assortment."),
    "9:3:3:1": ([9, 3, 3, 1], "Dihybrid cross with two traits showing independent segregation."),
    "12:3:1": ([12, 3, 1], "Dominant epistasis: dominant allele masks the effect of another."),
    "9:7": ([9, 7], "Duplicate recessive epistasis: both genes need dominant allele for trait."),
    "15:1": ([15, 1], "Duplicate dominant epistasis: one dominant allele from either gene is enough."),
    "9:3:4": ([9, 3, 4], "Recessive epistasis: recessive allele of one gene masks the other."),
    "13:3": ([13, 3], "Dominant & recessive (inhibitory) epistasis."),
    "9:6:1": ([9, 6, 1], "Polymeric gene interaction with additive effect.")
}

def test_segregation(observed, ratios):
    observed = np.array(observed)
    ratios = np.array(ratios)
    
    total = np.sum(observed)
    expected = ratios / np.sum(ratios) * total
    
    chi_stat, p_value = chisquare(f_obs=observed, f_exp=expected)
    
    return {
        'observed': observed.tolist(),
        'expected': expected.tolist(),
        'chi_square_stat': chi_stat,
        'p_value': p_value
    }

def compare_models(observed_counts):
    results = []
    for name, (ratio, explanation) in models.items():
        if len(ratio) == len(observed_counts):
            result = test_segregation(observed_counts, ratio)
            result['model'] = name
            result['ratio'] = ratio
            result['explanation'] = explanation
            results.append(result)
    if not results:
        return None
    best_result = max(results, key=lambda x: x['p_value'])
    return best_result

# UI layout
app_ui = ui.page_fluid(
    ui.h2("Genetic Segregation Ratio Tester"),
    ui.input_text_area("phenotypes", "Enter phenotypic values (comma-separated or column-paste)", rows=6),
    ui.output_ui("result_ui"),
    ui.output_plot("barplot")
)

# Server logic
def server(input, output, session):

    @reactive.Calc
    def observed_counts():
        text = input.phenotypes().strip()
        if not text:
            return None
        
        # Accept either comma-separated or newline-separated input
        raw_values = [v.strip() for v in text.replace(",", "\n").splitlines() if v.strip()]
        if not raw_values:
            return None
        df = pd.DataFrame({"phenotype": raw_values})
        counts = df["phenotype"].value_counts().sort_index()
        return counts

    @output
    @render.ui
    def result_ui():
        counts = observed_counts()
        if counts is None or len(counts) < 2:
            return ui.p("Please enter at least two phenotype classes.")

        result = compare_models(counts.values)
        if result is None:
            return ui.p("No segregation model matches the number of phenotype categories.")

        model_name = result['model']
        explanation = result['explanation']
        p_value = result['p_value']
        obs_list = result['observed']
        exp_list = [round(e, 2) for e in result['expected']]
        message = (
            f"The observed ratio {obs_list} does **not** significantly deviate from the expected ratio {exp_list} "
            f"(Chi-square p = {p_value:.4f}), so it is **consistent** with the best fit model: **{model_name}**."
            if p_value > 0.05 else
            f"The observed ratio {obs_list} **significantly deviates** from the expected ratio {exp_list} "
            f"(Chi-square p = {p_value:.4f}), so it is **not consistent** with the model: **{model_name}**."
        )

        return ui.card(
            ui.h4(f"Best Fitting Model: {model_name}"),
            ui.p(f"Observed counts: {obs_list}"),
            ui.p(f"Expected counts: {exp_list}"),
            ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
            ui.p(f"P-value: {p_value:.4f}"),
            ui.markdown(f"**Interpretation:** {message}"),
            ui.hr(),
            ui.markdown(f"**Model Explanation:** {explanation}")
        )

    @output
    @render.plot
    def barplot():
        counts = observed_counts()
        if counts is None:
