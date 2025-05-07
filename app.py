from shiny import App, ui, render, reactive
import pandas as pd
from collections import Counter
from scipy.stats import chisquare
import plotly.graph_objects as go

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

# Interpretation text for each model
MODEL_INTERPRETATIONS = {
    "3:1": "This pattern suggests a single gene with a dominant and recessive allele. One phenotype appears three times as often as the other, indicating classic Mendelian dominance.",
    "1:2:1": "This ratio indicates additive gene action in a monohybrid cross, where heterozygotes show an intermediate phenotype.",
    "12:3:1": "Dominant epistasis occurs when a dominant allele at one locus masks the expression of alleles at another locus, resulting in fewer phenotypic categories than expected.",
    "9:7": "This pattern is due to duplicate recessive epistasis, where both genes must have at least one dominant allele to produce a specific phenotype.",
    "15:1": "This is the result of duplicate dominant epistasis. A dominant allele at either of two loci produces the same phenotype, making the unique phenotype very rare.",
    "9:3:4": "Recessive epistasis involves one gene masking the expression of another, but only when it's homozygous recessive. This gives a characteristic 9:3:4 phenotypic distribution.",
    "13:3": "This ratio reflects a dominant and recessive interaction, where one dominant allele inhibits the expression of the other gene, leading to suppression of one phenotype.",
    "9:6:1": "This pattern, also known as polymeric gene interaction, shows additive effects of two genes producing a novel phenotype when both genes contribute equally.",
    "1:1:1:1": "This suggests independent assortment of two genes with complete dominance, typically seen in a dihybrid test cross.",
    "9:3:3:1": "A classic dihybrid Mendelian ratio involving two independently assorting genes with complete dominance."
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
        cleaned = [line.strip().lower() for line in lines]
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
        obs_counts = parsed_data()
        total = sum(obs_counts.values())
        sorted_categories = sorted(obs_counts.keys())

        results = []

        for name, ratio in matched_models().items():
            expected = sorted([total * r / sum(ratio) for r in ratio])
            observed = sorted([obs_counts.get(cat, 0) for cat in sorted_categories])

            if len(observed) != len(expected):
                continue
            chi2, p = chisquare(f_obs=observed, f_exp=expected)
            results.append((name, chi2, p, observed, expected))

        return sorted(results, key=lambda x: -x[2])  # sort by p-value desc

    @output
    @render.ui
    def result_ui():
        if not parsed_data():
            return ui.p("Paste data to begin.")

        if not matched_models():
            return ui.p("❌ No Mendelian model matches the number of phenotypic categories.")

        best = test_results()[0]
        name, chi2, p, obs, exp = best
        categories = sorted(parsed_data().keys())

        # Create bar plot
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=categories,
            y=obs,
            name="Observed",
            marker_color="lightgreen"
        ))
        fig.add_trace(go.Bar(
            x=categories,
            y=exp,
            name=f"Expected ({name})",
            marker_color="darkgreen"
        ))
        fig.update_layout(
            barmode="group",
            title="Observed vs Expected Counts",
            legend_title="Model Explanation",
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)"
        )

        return ui.div(
            ui.card(
                ui.div(
                    ui.h4("Best Fit Model", style="color: lightgreen;"),
                    ui.p(f"Model: **{name}**"),
                    ui.p(f"Observed counts: {obs}"),
                    ui.p(f"Expected counts: {[round(e, 2) for e in exp]}"),
                    ui.p(f"Chi-square statistic: {chi2:.4f}"),
                    ui.p(f"P-value: {p:.4f}"),
                ),
                style="box-shadow: 2px 2px 10px #ccc; padding: 1rem;"
            ),
            ui.output_plot("bar_plot"),
            ui.card(
                ui.div(
                    ui.h4("Model Interpretation", style="color: lightgreen;"),
                    ui.markdown(MODEL_INTERPRETATIONS.get(name, "No interpretation available for this model."))
                ),
                style="box-shadow: 2px 2px 10px #ccc; padding: 1rem; margin-top: 1rem;"
            )
        )

    @output
    @render.plot
    def bar_plot():
        if not test_results():
            return None
        name, _, _, obs, exp = test_results()[0]
        categories = sorted(parsed_data().keys())
        df = pd.DataFrame({
            "Category": categories,
            "Observed": obs,
            "Expected": exp
        })
        return df.plot.bar(x="Category", color=["lightgreen", "darkgreen"], rot=0).get_figure()

# Create and run the app
app = App(app_ui, server)
