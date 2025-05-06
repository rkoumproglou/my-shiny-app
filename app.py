from shiny import App, reactive, render, ui
import numpy as np
import plotly.graph_objects as go
from plotly.graph_objects import Figure
from scipy.stats import chisquare

# Define segregation models
models = {
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

def test_segregation(observed, ratios):
    observed = np.array(observed)
    ratios = np.array(ratios)
    
    total = np.sum(observed)
    expected = ratios / np.sum(ratios) * total
    
    chi_stat, p_value = chisquare(f_obs=observed, f_exp=expected)
    
    result = {
        'observed': observed.tolist(),
        'expected': expected.tolist(),
        'chi_square_stat': chi_stat,
        'p_value': p_value,
        'best_fit_ratio': ratios.tolist()
    }
    return result

def compare_models(observed_counts):
    results = []
    for name, ratio in models.items():
        if len(ratio) == len(observed_counts):
            result = test_segregation(observed_counts, ratio)
            result['model'] = name
            results.append(result)
    if not results:
        return None
    best_result = max(results, key=lambda x: x['p_value'])
    return best_result

# UI layout with CSS for drop shadow
app_ui = ui.page_fluid(
    ui.h2("Genetic Segregation Ratio Tester"),
    ui.input_text_area("counts", "Paste phenotypic values (one per line)", placeholder="e.g.\nred\nred\nred\norange\norange"),
    ui.output_ui("result_ui"),
    ui.tags.style("""
        .graph-container {
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        }
    """)
)

# Server logic
def server(input, output, session):

    @reactive.Calc
    def observed_counts():
        try:
            # Process the input as a list of phenotypic values
            phenotypes = [x.strip() for x in input.counts().splitlines() if x.strip()]
            # Count the occurrences of each phenotype
            unique_phenotypes = set(phenotypes)
            counts = [phenotypes.count(phenotype) for phenotype in unique_phenotypes]
            return counts, list(unique_phenotypes)
        except Exception:
            return None, None

    @output
    @render.ui
    def result_ui():
        counts, unique_phenotypes = observed_counts()
        if not counts:
            return ui.p("Please enter valid phenotypic values (one per line).")
        
        result = compare_models(counts)
        if result is None:
            return ui.p("No model matches the length of the observed data.")
        
        # Determine p-value interpretation
        if result['p_value'] < 0.05:
            interpretation = f"The observed ratio does not significantly deviate from the expected ratio based on the Chi-square test, so it is considered consistent with the best fit model {result['model']}*."
        else:
            interpretation = "No model explains the observed segregation ratio based on the Chi-square test."
        
        # Construct the output
        output = ui.panel_well(
            ui.h4("Best-Fit Model Report", style="background-color: lightgreen; padding: 10px;"),
            ui.p(interpretation),
            ui.p(f"Best Fit Model: {result['model']}*"),
            ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
            ui.p(f"P-value: {result['p_value']:.4f}"),
            ui.p(f"Tested Ratio: {result['best_fit_ratio']}")
        )

        # Add graph if p-value is greater than 0.05
        if result['p_value'] >= 0.05:
            output += ui.card(
                ui.plotlyOutput("graph", height="400px"), 
                style="box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);"
            )

            # Explanation of each model for the best fit model
            output += ui.p(f"* {result['model']}* Model Explanation:"),
            output += ui.p(
                "* **3:1 Ratio**: This occurs when a single gene has one dominant and one recessive allele. F₂ Result: 3 = dominant phenotype, 1 = recessive phenotype. Genotypic breakdown: 1 homozygous dominant, 2 heterozygous, 1 homozygous recessive."
            )
            output += ui.p(
                "* **1:2:1 Ratio**: This ratio arises when there is incomplete dominance or codominance. F₂ Result: 1 = homozygous dominant phenotype, 2 = heterozygous intermediate or combined phenotype, 1 = homozygous recessive phenotype."
            )
            output += ui.p(
                "* **Recessive Epistasis (9:3:4)**: The recessive allele of one gene masks the expression of the other gene. 9 = both dominant alleles present → first phenotype, 3 = one dominant allele from the second gene only → second phenotype, 4 = homozygous recessive in the epistatic gene → third phenotype."
            )
            output += ui.p(
                "* **Dominant Epistasis (12:3:1)**: The dominant allele of one gene masks the expression of the other gene. 12 = presence of dominant epistatic allele → first phenotype, 3 = non-epistatic gene expressed in absence of dominant epistasis → second phenotype, 1 = both genes homozygous recessive → third phenotype."
            )
            output += ui.p(
                "* **Dominant and Recessive (Inhibitory) Epistasis (13:3)**: A phenotype appears either when a dominant inhibitor is present or when the other gene is homozygous recessive. 13 = phenotype suppressed due to dominant or recessive epistasis, 3 = phenotype expressed only when inhibitor is absent and the second gene has at least one dominant allele."
            )
            output += ui.p(
                "* **Duplicate Recessive Epistasis (9:7)**: Both genes must carry at least one dominant allele for the phenotype to be expressed. 9 = phenotype expressed, 7 = any homozygous recessive condition in either gene blocks expression."
            )
            output += ui.p(
                "* **Duplicate Dominant Epistasis (15:1)**: A dominant allele in either gene is sufficient to produce the phenotype. 15 = phenotype expressed when at least one dominant allele is present in either gene, 1 = alternate phenotype expressed only in the double recessive condition."
            )
            output += ui.p(
                "* **Polymeric Gene Interaction (9:6:1)**: Both genes contribute additively to the phenotype, and their combined presence produces an enhanced or new trait. 9 = both dominant alleles present → enhanced phenotype, 6 = only one dominant allele present → intermediate phenotype, 1 = both genes recessive → basic phenotype."
            )

        return output

    @output
    @render.plotly
    def graph():
        counts, unique_phenotypes = observed_counts()
        if not counts:
            return go.Figure()

        result = compare_models(counts)
        if result is None or result['p_value'] >= 0.05:
            return go.Figure(
                data=[
                    go.Bar(
                        x=unique_phenotypes,
                        y=result['observed'],
                        name='Observed',
                        marker=dict(color="rgb(204,204,255)")
                    ),
                    go.Bar(
                        x=unique_phenotypes,
                        y=result['expected'],
                        name='Expected',
                        marker=dict(color="rgb(255,204,204)")
                    ),
                ],
                layout=go.Layout(
                    title="Observed vs. Expected Distribution",
                    showlegend=True,
                    xaxis={"title": "Phenotypes"},
                    yaxis={"title": "Counts"},
                    bargap=0.3  # To prevent the bars from touching each other
                )
            )
        return go.Figure()

# Create app
app = App(app_ui, server)
