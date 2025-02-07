using CSV
using DataFrames
using QuartetNetworkGoodnessFit
using PhyloNetworks
using RCall
using Distributions
using GenieFramework
using PlotlyBase
@genietools

using Genie, Genie.Router, Genie.Renderer.Html, Genie.Requests

const species_tree = readTopology("(Smi165:18.12505298,(Age001:8.747492442,(Adi001:8.308140028,(((Asu001:3.171531094,Aga001:3.171531094):2.767496398,(Ape001:1.381629136,(Aza037:1.363263922,Ama006:1.363263922):0.01836521346):4.555985647):0.1963665132,Aru001:6.135394005):2.172746023):0.4393524143):9.377560535);")
const gene = joinpath("uploads", "gene_tree.rtf")
genetrees = readMultiTopology(gene)
q,t = countquartetsintrees(genetrees)
df = writeTableCF(q,t)
CSV.write("tableCF.csv", df)
const geneCF = readTableCF("tableCF.csv")
topologyMaxQPseudolik!(species_tree, geneCF);
df_long = fittedQuartetCF(geneCF, :long)

# parameter values
const TFchoice = [true, false]
const stat_choice = [Symbol("LRT"), Symbol("Qlog"), Symbol("pearson")]
const corr_choice = [Symbol("simulation"), Symbol("none")]

mean_Norm = 0
stddev = 1

x1_range = range(-3*stddev, 3*stddev, length=100)
y1 = pdf(Normal(mean_Norm, stddev), x1_range)

res = quarnetGoFtest!(species_tree, geneCF,true;quartetstat=:LRT,correction=:simulation,seed=100,nsim=100,verbose=false,keepfiles=false)
overall_p = res[1]
uncorrected_z = res[2]
sigma = res[3]
pval_vec = res[4]
optnet = res[5]
z_vec = res[6]
x2_range = range(-3*sigma, 3*sigma, length=100)
y2 = pdf(Normal(mean_Norm, sigma), x2_range)

@handlers begin
    @out TFchoice
    @out stat_choice
    @out corr_choice

    @in optimal_bl = true
    @in test_stat = Symbol("LRT")
    @in correct = Symbol("simulation")
    @in setseed = 100
    @in nsimulation = 100
    @in diagnose = false
    @in keep = false
    
    @out Pval = overall_p
    @out trace_cf = [PlotlyBase.scatter(x=df_long[:, :obsCF],y=df_long[:, :expCF], mode = "markers")]
    @out layout_cf = PlotlyBase.Layout(title="CF Scatter Plot", 
    xaxis=attr(title="Observed CF"), yaxis=attr(title="Expected CF"))
    @out trace_p = [histogram(x=pval_vec, name="P Value Distribution")]
    @out layout_p = PlotlyBase.Layout(
        xaxis_title="P-value",
        yaxis_title="Count")
    @out Z_combined_trace = [histogram(x=z_vec, name="Z-Value Distribution", histnorm="probability density", opacity=0.5),
    scatter(x=x1_range, y=y1, mode="lines", name="Standard Normal", line=attr(color="red", width=2)),
    scatter(x=x2_range, y=y2, mode="lines", name="Fitted Normal", line=attr(color="blue", width=2))]
    @out Zhist_layout = PlotlyBase.Layout(
        title="Z-Value Distribution",
        xaxis_title="Z-value",
        yaxis_title="Density",
        barmode="overlay"
    )

    @in isready = false
    @onchange isready,optimal_bl,test_stat,correct,setseed,nsimulation,diagnose,keep begin
        res = quarnetGoFtest!(species_tree, geneCF,optimal_bl;quartetstat=test_stat,correction=correct,seed=setseed,nsim=nsimulation,verbose=diagnose,keepfiles=keep)
        overall_p = res[1]
        uncorrected_z = res[2]
        sigma = res[3]
        pval_vec = res[4]
        optnet = res[5]
        z_vec = res[6]
        x2_range = range(-3*sigma, 3*sigma, length=100)
        y2 = pdf(Normal(mean_Norm, sigma), x2_range)
        
        Pval = overall_p
        trace_p = [histogram(x=pval_vec, name="P Value Distribution")]
        Z_combined_trace = [histogram(x=z_vec, name="Z-Value Distribution", histnorm="probability density", opacity=0.5),
            scatter(x=x1_range, y=y1, mode="lines", name="Standard Normal", line=attr(color="red", width=2)),
            scatter(x=x2_range, y=y2, mode="lines", name="Fitted Normal", line=attr(color="blue", width=2))]

        layout_p = PlotlyBase.Layout(
                xaxis_title="P-value",
                yaxis_title="Count")
        Zhist_layout = PlotlyBase.Layout(
            title="Z-Value Distribution",
            xaxis_title="Z-value",
            yaxis_title="Density",
            barmode="overlay"
        )
    end
end

@page("/","ui.html")