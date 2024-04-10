# expect CF
# obs CF
# scatter plot

using CSV
using DataFrames
using QuartetNetworkGoodnessFit
using PhyloNetworks
using Plots

## reading files from users
gene_CF = DataFrame(CSV.File(joinpath(dirname(pathof(QuartetNetworkGoodnessFit)), "..","test","example_qCF_5taxa.csv")), copycols=false);
## large data for test:
# DataFrame(CSV.File(joinpath(dirname(pathof(QuartetNetworkGoodnessFit)), "..","test","larger_test_data.csv")), copycols=false);
## smaller data for test:
# DataFrame(CSV.File(joinpath(dirname(pathof(QuartetNetworkGoodnessFit)), "..","test","example_qCF_5taxa.csv")), copycols=false); # gene structures - later input from users
tree = readTopology("((((D,C),(A,B)),E),O);")
## large tree for test:
# readTopology("(Smic:18.12505298,(Agre:8.747492442,(Adig:8.308140028,(((Asua:3.171531094,Agra:3.171531094):2.767496398,(Aper:1.381629136,(Azaa:1.363263922,Amad:1.363263922):0.01836521346):4.555985647):0.1963665132,Arub:6.135394005):2.172746023):0.4393524143):9.377560535);")
## smaller tree for test:
# readTopology("((((D,C),(A,B)),E),O);") # tree structure - input from users
# this tree struture now is the correct one


## set parameters
# need to make choices for users to choose
optimal_bl = true   # optimize the branch length or not
test_stat = :LRT
correct = :simulation   # correct the dependence of 4-taxon
set_seed = 2
num_sim = 1000  # number of simulations
diagnose = false    # if we want to see the simulation process
keep = false    # choose to keep the simulation file or not


## run the function
res = quarnetGoFtest!(tree, gene_CF, optimal_bl;)

## see the simulation process
res_see_sim = quarnetGoFtest!(tree, gene_CF, optimal_bl; verbose = true)

## result without correction
res_no_corr = quarnetGoFtest!(tree, gene_CF, optimal_bl; correction = :none)
# not be able to get corrected z vectors and no simulation.

# ## Meaning of the contents in res
#
# It contains 6 outputs:
# - overall p-value: a number
# - uncorrected z: a number
# - $\sigma$: a number
# - a vector of p-values for each 4-taxon
# - the network with optimal branch length
# - a vector of empirical z values

# +
# separate the output res
overall_p = res[1]
uncorrected_z = res[2]
sigma = res[3]
pval_vec = res[4]
optnet = res[5]
z_vec = res[6]

length(z_vec) # consistent with the number of simulations
# compute the proportion of empirical z values that are greater
# then the uncorrected_z
count = 0
for z in z_vec
    if z > uncorrected_z
        count += 1
    elseif z == uncorrected_z
        count += 1
    end
end

# print(count)
# divide count from the length of z_vec to compute the proportion
prop_greater = count / length(z_vec)
print("The overall p-value is: ", overall_p, "\n")
print("The empirical p-value is: ", prop_greater, "\n")
print("The overall p-value rejects the species tree.", "\n")
print("The empirical p-value does not reject?")
# -

# # visualize the result
# we want to visualize the p-value, z-value, sigma
# maybe later we want to know which are the genes that reject the tree, and give an alignmnet.

print(pval_vec)

# +
## plot the p-value vector for each 4-taxon
# create a dataframe of index of 4-taxon and p-value
# index as x, p-value as y axis
# idx = 1:15
# bar(idx, pval_vec)

# range(0, 1, length = 7)
histogram(pval_vec, bins = range(0, 1, length = 21))  # choose 21 to have 20 bars b/c 20 will cut 0.05 out.
## add the 0.05 cutoff
vline!([0.05], lw = 3)

# +
## let's see the proportional plot of p-values

# # create a new dataframe with count of p-values > 0.05 and < 0.05
# pval_cate = ["< 0.05", ">= 0.05"]

# #=
# write a function to count the p-values
# it takes the p-value vector as input
# counts the p-values that are < 0.05, and >= 0.05
# return a count vector
# =#
# function count_pval(p_vec)
#     count = [0,0]
#     for p in p_vec
#         if p < 0.05
#             count[1] += 1
#         elseif p >= 0.05
#             count[2] += 1
#         end
#     end
#     return count
# end

# # define count vector
# count_vec = count_pval(pval_vec)



# +
# combine count vector and the category of pvalue into a df
# prop_df = DataFrame("pval_category" => pval_cate, "count" => count_vec)

# +
# plot(prop_df)
# -

# ### Z value

# shape of z-value vector
size(z_vec) # 1000

# We have 1000 simulations, each simulation contains info that has the same species tree and same number of gene trees. In each simulation, repeat the procedure to find the outlier p-values, and the proportion of p-values < 0.05. This prop is the $\hat{p}_{out}$ in the formula to compute z. So we can obtain $z_i$ in each simulation.

# +
# bar(z_vec)
# -

# **Want:** z-value and probability (frequency)

# +
unique_value = unique(z_vec) # a vector with unique z-values
count_z = zeros(length(unique_value))
for j in 1:length(unique_value)
    for i in 1:length(z_vec)
        if z_vec[i] == unique_value[j]
            count_z[j] += 1
        end
    end
end

print(count_z)

# +
# use the count vector to compute the proportion of each z
# length(z_vec) is the total number of elements in z_vec
prop_z = zeros(length(count_z))
for i in 1:length(count_z)
    prop_z[i] = count_z[i]/length(z_vec)
end

print(prop_z)

# +
using Distributions
d = Normal(0,1)

bar(prop_z)
plot!(x -> pdf(d, x), lw = 3)
# -

# This is a plot showing the distribution of z and compare with the standard normal.

# What if we want to see the comparison with the normal with new sigma?

# +
new_normal = Normal(0, sigma)
# Normal function takes mu and sigma

bar(prop_z)
plot!(x -> pdf(new_normal, x), lw = 3)
# -

print(unique_value)

histogram(z_vec, normalize = :pdf)
plot!(x -> pdf(d, x), lw = 3, color = :red)
plot!(x -> pdf(new_normal, x), lw = 3, color = :green)

# Can we say that the new normal distribution is a better approximation of z random variables from the plot?
#
# I think: yes, but this data may not be a good choice to see this, so we want to see the plots with a larger data.

# ### Uncorrected z-value

# ## Compare the optimized tree
#
# Can we make it prettier?

using PhyloPlots

# plot the original species tree
PhyloPlots.plot(tree, useedgelength = true)

# plot the optimal branch tree
PhyloPlots.plot(optnet, useedgelength = true, showedgelength = true)

# ## Change test stat

res1 = quarnetGoFtest!(tree, gene_CF, optimal_bl; quartetstat=:Qlog)

overall_p1 = res1[1]
uncorrected_z1 = res1[2]
sigma1 = res1[3]
pval_vec1 = res1[4]
optnet1 = res1[5]
z_vec1 = res1[6]

overall_p1

# ### p-value plot for Qlog

histogram(pval_vec1, bins = range(0, 1, length = 21))  # choose 21 to have 20 bars b/c 20 will cut 0.05 out.
## add the 0.05 cutoff
vline!([0.05], lw = 3)

print(pval_vec)

print(pval_vec1)

# ### z-value distribution for Qlog

print(sigma, "\n")
print(sigma1)

# +
new_normal1 = Normal(0, sigma1)

histogram(z_vec1, normalize = :pdf)
plot!(x -> pdf(d, x), lw = 3, color = :red)
plot!(x -> pdf(new_normal1, x), lw = 3, color = :green)
# -

# ### trees
#
# No change after changing the test stat

PhyloPlots.plot(optnet1)

# ## pearson statistics

res2 = quarnetGoFtest!(tree, gene_CF, optimal_bl; quartetstat=:pearson)

overall_p2 = res2[1]
uncorrected_z2 = res2[2]
sigma2 = res2[3]
pval_vec2 = res2[4]
optnet2 = res2[5]
z_vec2 = res2[6]

overall_p2

# ### p-value for pearson

histogram(pval_vec2, bins = range(0, 1, length = 21))  # choose 21 to have 20 bars b/c 20 will cut 0.05 out.
## add the 0.05 cutoff
vline!([0.05], lw = 3)

# ### z-value distribution for pearson

# +
new_normal2 = Normal(0, sigma2)

histogram(z_vec2, normalize = :pdf)
plot!(x -> pdf(d, x), lw = 3, color = :red)
plot!(x -> pdf(new_normal2, x), lw = 3, color = :green)
