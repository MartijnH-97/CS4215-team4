import numpy as np
from scipy.stats import f, t
import math
from LaTeXPrinter import create_table
from dataCollector import collector
import matplotlib.pyplot as plt

# Toggle creation of output graphs.
create_output = True

# Define confidence level to use.
alpha = 0.10


def square(x): return x ** 2


# Obtain data information
titles, names, parameters, DATA = collector()
a, b, c, r = parameters

# Calculate averages
AVERAGES = np.zeros((a, c, b))

for i in range(0, a):
    for j in range(0, b):
        for k in range(0, c):
            average = np.average(DATA[i][k*r : (k+1)*r], axis=0)
            AVERAGES[i][k] = average

SUMS = [[], [], []]
MEANS = [[], [], []]

# Calculate sum and mean of each level of factor A
for level in range(0, a):
    SUMS[0].append(np.sum(AVERAGES[level]))
    MEANS[0].append(np.average(AVERAGES[level]))

# Calculate sum and mean of each level of factor B
for level in range(0, b):
    partialSum = 0.0
    for depth in range(0, a):
        partialSum += np.sum(AVERAGES[depth], 0)[level]
    SUMS[1].append(partialSum)
    MEANS[1].append(partialSum/(a*c))

# Calculate sum and mean of each level of factor C
for level in range(0, c):
    partialSum = 0.0
    for depth in range(0, a):
        partialSum += np.sum(AVERAGES[depth], 1)[level]
    SUMS[2].append(partialSum)
    MEANS[2].append(partialSum/(a*b))

# # Calculate total sum and mean
TOTAL_SUM = np.sum(AVERAGES)
TOTAL_MEAN = np.average(AVERAGES)

# Calculate means of different combinations
MEANS_AB = np.zeros((a, b))
for i in range(0, a):
    for j in range(0, b):
        sum = 0.0
        for k in range(0, c):
            sum += AVERAGES[i][k][j]

        MEANS_AB[i][j] = sum / c

MEANS_AC = np.zeros((a, c))
for i in range(0, a):
    for k in range(0, c):
        sum = 0.0
        for j in range(0, b):
            sum += AVERAGES[i][k][j]

        MEANS_AC[i][k] = sum / b

MEANS_BC = np.zeros((b, c))
for j in range(0, b):
    for k in range(0, c):
        sum = 0.0
        for i in range(0, a):
            sum += AVERAGES[i][k][j]

        MEANS_BC[j][k] = sum / a

# Calculate main effects
EFFECTS = map(lambda x: x - TOTAL_MEAN, MEANS)

# Calculate the 2-way interactions
INTERACTIONS_AB = np.zeros((a, b))
for i in range(0, a):
    for j in range(0, b):
        INTERACTIONS_AB[i][j] = MEANS_AB[i][j] - MEANS[0][i] - MEANS[1][j] + TOTAL_MEAN

INTERACTIONS_AC = np.zeros((a, c))
for i in range(0, a):
    for k in range(0, c):
        INTERACTIONS_AC[i][k] = MEANS_AC[i][k] - MEANS[0][i] - MEANS[2][k] + TOTAL_MEAN

INTERACTIONS_BC = np.zeros((b, c))
for j in range(0, b):
    for k in range(0, c):
        INTERACTIONS_BC[j][k] = MEANS_BC[j][k] - MEANS[1][j] - MEANS[2][k] + TOTAL_MEAN

INTERACTIONS_ABC = np.zeros((a, c, b))
for k in range(0, c):
    for j in range(0, b):
        for i in range(0, a):
            INTERACTIONS_ABC[i][k][j] = AVERAGES[i][k][j] - MEANS_AB[i][j] - MEANS_AC[i][k] - MEANS_BC[j][k] + MEANS[0][i] + MEANS[1][j] + MEANS[2][k] - TOTAL_MEAN

# Create store for SS values
SS = np.zeros((8, 1))
# Calculate Sum of Squares for A, B and C
SS[0] = c*b*r*np.sum(square(EFFECTS[0]))
SS[1] = c*a*r*np.sum(square(EFFECTS[1]))
SS[2] = a*b*r*np.sum(square(EFFECTS[2]))

# Calculate Sum of Squares for 2-way interactions
sums = 0.0
for i in range(0, a):
    for j in range(0, b):
        sums += INTERACTIONS_AB[i][j]**2
SS[3] = r*c*sums

sums = 0.0
for i in range(0, a):
    for k in range(0, c):
        sums += INTERACTIONS_AC[i][k]**2
SS[4] = r*b*sums

sums = 0.0
for j in range(0, b):
    for k in range(0, c):
        sums += INTERACTIONS_BC[j][k]**2
SS[5] = r*a*sums

sums = 0.0
for i in range(0, a):
    for j in range(0, b):
        for k in range(0, c):
            sums += INTERACTIONS_ABC[i][k][j]**2
SS[6] = r*sums

SSY = np.sum(np.sum(np.sum(square(DATA))))
SSO = a * b * c * r * (TOTAL_MEAN ** 2)
SST = SSY - SSO
SS[7] = SST - np.sum(SS)

# Store percentage of variation explained
percentage_variation_exp = np.zeros((8, 1))
percentage_variation_exp[0] = SS[0]/SST                 # A
percentage_variation_exp[1] = SS[1]/SST                 # B
percentage_variation_exp[2] = SS[2]/SST                 # C
percentage_variation_exp[3] = SS[3]/SST                 # AB
percentage_variation_exp[4] = SS[4]/SST                 # AC
percentage_variation_exp[5] = SS[5]/SST                 # BC
percentage_variation_exp[6] = SS[6]/SST                 # ABC
percentage_variation_exp[7] = SS[7]/SST                 # Error

# Transform from fractions to percents
percentage_variation_exp = percentage_variation_exp*100

# Store degrees of freedom
dof = np.zeros((8, 1))
dof[0] = a - 1                  # A
dof[1] = b - 1                  # B
dof[2] = c - 1                  # C
dof[3] = dof[0]*dof[1]          # AB
dof[4] = dof[0]*dof[2]          # AC
dof[5] = dof[1]*dof[2]          # BC
dof[6] = dof[0]*dof[1]*dof[2]   # ABC
dof[7] = a*b*c*(r - 1)          # Error

# Store MS values
MS = np.zeros((8, 1))

MS[0] = SS[0]/dof[0]
MS[1] = SS[1]/dof[1]
MS[2] = SS[2]/dof[2]
MS[3] = SS[3]/dof[3]
MS[4] = SS[4]/dof[4]
MS[5] = SS[5]/dof[5]
MS[6] = SS[6]/dof[6]
MS[7] = SS[7]/dof[7]

# Calculate F-comp values
Fcomp = np.zeros((7, 1))
Fcomp[0] = MS[0]/MS[7]
Fcomp[1] = MS[1]/MS[7]
Fcomp[2] = MS[2]/MS[7]
Fcomp[3] = MS[3]/MS[7]
Fcomp[4] = MS[4]/MS[7]
Fcomp[5] = MS[5]/MS[7]
Fcomp[6] = MS[6]/MS[7]

# Calculate F-table values
Ftab = np.zeros((7, 1))
Ftab[0] = f.ppf(1 - alpha, dof[0], dof[7], loc=0, scale=1)
Ftab[1] = f.ppf(1 - alpha, dof[1], dof[7], loc=0, scale=1)
Ftab[2] = f.ppf(1 - alpha, dof[2], dof[7], loc=0, scale=1)
Ftab[3] = f.ppf(1 - alpha, dof[3], dof[7], loc=0, scale=1)
Ftab[4] = f.ppf(1 - alpha, dof[4], dof[7], loc=0, scale=1)
Ftab[5] = f.ppf(1 - alpha, dof[5], dof[7], loc=0, scale=1)
Ftab[6] = f.ppf(1 - alpha, dof[6], dof[7], loc=0, scale=1)
#
# print Ftab

s_e = math.sqrt(MS[7])

standardDeviations = []

standardDevA = []
for i in range(0, a):
    standardDevA.append(s_e*math.sqrt(dof[0]/(a*b*c*r)))
standardDeviations.append(standardDevA)

standardDevB = []
for j in range(0, b):
    standardDevB.append(s_e*math.sqrt(dof[1]/(a*b*c*r)))
standardDeviations.append(standardDevB)

standardDevC = []
for k in range(0, c):
    standardDevC.append(s_e*math.sqrt(dof[2]/(a*b*c*r)))
standardDeviations.append(standardDevC)

# Confidence intervals for main effects
CIs_MAIN = []

CI_A = []
for i in range(0, a):
    high = EFFECTS[0][i] + t.ppf(1 - alpha, dof[7])*standardDeviations[0][i]
    low = EFFECTS[0][i] - t.ppf(1 - alpha, dof[7])*standardDeviations[0][i]
    CI_A.append((low, high))
CIs_MAIN.append(CI_A)

CI_B = []
for j in range(0, b):
    high = EFFECTS[1][j] + t.ppf(1 - alpha, dof[7])*standardDeviations[1][j]
    low = EFFECTS[1][j] - t.ppf(1 - alpha, dof[7])*standardDeviations[1][j]
    CI_B.append((low, high))
CIs_MAIN.append(CI_B)

CI_C = []
for k in range(0, c):
    high = EFFECTS[2][k] + t.ppf(1 - alpha, dof[7])*standardDeviations[2][k]
    low = EFFECTS[2][k] - t.ppf(1 - alpha, dof[7])*standardDeviations[2][k]
    CI_C.append((low, high))
CIs_MAIN.append(CI_C)

# Due to being a 2^3 design, all standard deviations are the same.
STD = standardDeviations[0][0]

# OUTPUT GENERATION STARTS HERE

# #  Create LaTeX table of the ANOVA results.
RESULTS = np.zeros((11, 6))

# Set sum of squares field
RESULTS[0:3, 0] = [SSY, SSO, SST]
RESULTS[3:, 0] = SS[:, 0]

# Set Percentage of variation field
RESULTS[2, 1] = 100
RESULTS[3:, 1] = percentage_variation_exp[:, 0]

# Set Degree of Freedom field
RESULTS[3:, 2] = dof[:, 0]

# Set Mean Square field
RESULTS[3:, 3] = MS[:, 0]

# Set F-Comp field
RESULTS[3:-1, 4] = Fcomp[:, 0]

# Set F-Table field
RESULTS[3:-1, 5] = Ftab[:, 0]


table_row_names = [
                    "$y$",
                    "$\\bar{y}...$",
                    "$y - \\bar{y}...$"] + \
                    [titles[0]] + \
                    [titles[1]] + \
                    [titles[2]] + \
                    [titles[0] + " x "+ titles[1]] + \
                    [titles[0] + " x "+ titles[2]] + \
                    [titles[1] + " x "+ titles[2]] + \
                    [titles[0] + " x "+ titles[1] + " x " + titles[2]] + \
                    ["error"]
table_col_names = ["Component", "Sum of squares", "\\% Variation", "DF", "Mean Square", "F-Comp.", "F-Table"]

title = "ANOVA results"
label = "anova"
create_table(RESULTS, table_row_names,  table_col_names, title, label)

# #  Create LaTeX table of the averages used during the study.
average_table_data = np.zeros((a * c + a, b))
average_table_data[1, :] = AVERAGES[0, 0, :]
average_table_data[2, :] = AVERAGES[0, 1, :]
average_table_data[4, :] = AVERAGES[1, 0, :]
average_table_data[5, :] = AVERAGES[1, 1, :]

table_row_names = [
                        titles[0] + " = " + str(names[0][0]),
                        titles[2] + " = " + str(names[2][0]),
                        titles[2] + " = " + str(names[2][1]),
                        titles[0] + " = " + str(names[0][1]),
                        titles[2] + " = " + str(names[2][0]),
                        titles[2] + " = " + str(names[2][1])
]
table_col_names = ["Factors", titles[1] + " = " + str(names[1][0]), titles[1] + " = " + str(names[1][1])]

title = "Averages over the repetitions of the collected data"
label = "averages"
create_table(average_table_data, table_row_names, table_col_names, title, label)

# #  Create LaTeX table of the raw data.
all_data_table_data = np.zeros((a * c * r + a, b))
all_data_table_data[1: 1+a * r] = DATA[0, :]
all_data_table_data[2 + a * r : 2 + 2*a * r] = DATA[1, :]

table_row_names = [
    titles[0] + " = " + str(names[0][0]),
    titles[2] + " = " + str(names[2][0]),
    "",
    "",
    titles[2] + " = " + str(names[2][1]),
    "",
    "",
    titles[0] + " = " + str(names[0][1]),
    titles[2] + " = " + str(names[2][0]),
    "",
    "",
    titles[2] + " = " + str(names[2][1]),
    "",
    "",
]

table_col_names = ["Factors", titles[1] + " = " + str(names[1][0]), titles[1] + " = " + str(names[1][1])]

title = "Data collected during the study"
label = "rawdata"
create_table(all_data_table_data, table_row_names, table_col_names, title, label)

# Create LaTex table for the confidence intervals of the main effects.
STD_mu = s_e*math.sqrt(1.0 / (a*b*r*c))
CI_mu = [TOTAL_MEAN + t.ppf(1 - alpha, dof[7])*STD_mu, TOTAL_MEAN - t.ppf(1 - alpha, dof[7])*STD_mu]

effects_table_data = [[TOTAL_MEAN, STD_mu, "({} , {})".format(CI_mu[0][0], CI_mu[1][0])]]

effects_table_data = np.append(effects_table_data, [["", "", ""]], axis=0)
for i in range(0, len(CIs_MAIN[0])):
    effects_table_data = np.append(effects_table_data, [[EFFECTS[0][i], standardDeviations[0][i], "({} , {})".format(CIs_MAIN[0][i][0][0], CIs_MAIN[0][i][1][0])]], axis=0)

effects_table_data = np.append(effects_table_data, [["", "", ""]], axis=0)
for j in range(0, len(CIs_MAIN[1])):
    effects_table_data = np.append(effects_table_data, [[EFFECTS[1][j], standardDeviations[1][j], "({} , {})".format(CIs_MAIN[1][j][0][0], CIs_MAIN[1][j][1][0])]], axis=0)

effects_table_data = np.append(effects_table_data, [["", "", ""]], axis=0)
for k in range(0, len(CIs_MAIN[2])):
    effects_table_data = np.append(effects_table_data, [[EFFECTS[2][k], standardDeviations[2][k], "({} , {})".format(CIs_MAIN[2][k][0][0], CIs_MAIN[2][k][1][0])]], axis=0)

table_row_names = ["$\\mu$"] + [titles[0]] + names[0] + [titles[1]] + names[1] + [titles[2]] + names[2]
table_col_names = ["Parameter", "Mean Effect", "Standard Deviation", "CI"]

title = "Confidence intervals of main effects"
label = "CI_effects"
create_table(effects_table_data, table_row_names, table_col_names, title, label)

# Create LaTeX table of the AB interaction.
g=0
h=1
effectsAB_table_data = [[titles[h] + " = " + str(names[h][0]), "", "", titles[h] + " = " + str(names[h][1]), "", ""]]
for i in range(0, a):
    frst_effect = INTERACTIONS_AB[i][0]
    bar = t.ppf(1 - alpha, dof[7])*STD

    frst_CI = "({} ,{})".format(frst_effect + bar[0], frst_effect - bar[0])

    scnd_effect = INTERACTIONS_AB[i][1]
    scnd_CI = "({} ,{})".format(scnd_effect + bar[0], scnd_effect - bar[0])

    row = [[frst_effect, STD, frst_CI, scnd_effect, STD, scnd_CI]]
    effectsAB_table_data = np.append(effectsAB_table_data, row, axis=0)

table_row_names = ["", titles[g] + " = " + str(names[g][0])] + [titles[g] + " = " + str(names[g][1])]
table_col_names = ["Parameter",
                   "Mean Effect interaction",
                   "Standard Deviation",
                   "CI",
                   "Mean Effect interaction",
                   "Standard Deviation",
                   "CI"]

title = "Effects and confidence intervals of interaction between " + titles[g] + " and " + titles[h]
label = "CI_effects_interactions_AB"
create_table(effectsAB_table_data, table_row_names, table_col_names, title, label)

# Create LaTeX table of the AC interaction.
g=0
h=2
effectsAC_table_data = [[titles[h] + " = " + str(names[h][0]), "", "", titles[h] + " = " + str(names[h][1]), "", ""]]
for i in range(0, a):
    frst_effect = INTERACTIONS_AC[i][0]
    bar = t.ppf(1 - alpha, dof[7])*STD

    frst_CI = "({} ,{})".format(frst_effect + bar[0], frst_effect - bar[0])

    scnd_effect = INTERACTIONS_AC[i][1]
    scnd_CI = "({} ,{})".format(scnd_effect + bar[0], scnd_effect - bar[0])

    row = [[frst_effect, STD, frst_CI, scnd_effect, STD, scnd_CI]]
    effectsAC_table_data = np.append(effectsAC_table_data, row, axis=0)

table_row_names = ["", titles[g] + " = " + str(names[g][0])] + [titles[g] + " = " + str(names[g][1])]
table_col_names = ["Parameter",
                   "Mean Effect interaction",
                   "Standard Deviation",
                   "CI",
                   "Mean Effect interaction",
                   "Standard Deviation",
                   "CI"]

title = "Effects and confidence intervals of interaction between " + titles[g] + " and " + titles[h]
label = "CI_effects_interactions_AC"
create_table(effectsAC_table_data, table_row_names, table_col_names, title, label)

# Create LaTeX table of the BC interaction.
g=1
h=2
effectsAB_table_data = [[titles[h] + " = " + str(names[h][0]), "", "", titles[h] + " = " + str(names[h][1]), "", ""]]
for i in range(0, a):
    frst_effect = INTERACTIONS_BC[i][0]
    bar = t.ppf(1 - alpha, dof[7])*STD

    frst_CI = "({} ,{})".format(frst_effect + bar[0], frst_effect - bar[0])

    scnd_effect = INTERACTIONS_BC[i][1]
    scnd_CI = "({} ,{})".format(scnd_effect + bar[0], scnd_effect - bar[0])

    row = [[frst_effect, STD, frst_CI, scnd_effect, STD, scnd_CI]]
    effectsAB_table_data = np.append(effectsAB_table_data, row, axis=0)

table_row_names = ["", titles[g] + " = " + str(names[g][0])] + [titles[g] + " = " + str(names[g][1])]
table_col_names = ["Parameter",
                   "Mean Effect interaction",
                   "Standard Deviation",
                   "CI",
                   "Mean Effect interaction",
                   "Standard Deviation",
                   "CI"]

title = "Effects and confidence intervals of interaction between " + titles[g] + " and " + titles[h]
label = "CI_effects_interactions_BC"
create_table(effectsAB_table_data, table_row_names, table_col_names, title, label)

# Create plots of fitted value against residual and fitted value against actual value.
fitted_residuals = []
fitted_actual = []
for i in range(0, a):
    for k in range(0, c):
        block = DATA[i][k*r:(k+1)*r]
        for block_row in block:
            for j in range(0, b):
                actual = block_row[j]
                main_effects = EFFECTS[0][i] + EFFECTS[1][j] + EFFECTS[2][k]
                second_effects = INTERACTIONS_AB[i][j] + INTERACTIONS_AC[i][k] + INTERACTIONS_BC[j][k]
                fitted = TOTAL_MEAN + main_effects + second_effects + INTERACTIONS_ABC[i][k][j]
                residual = actual - fitted

                fitted_residuals.append((fitted, residual))
                fitted_actual.append((fitted, actual))

title = "Residual plots"
plt.figure()
plt.title(title)
plt.subplot(211)
plt.scatter(*zip(*fitted_residuals))
plt.xlabel('Predicted value')
plt.ylabel('Residual')

plt.subplot(212)
plt.scatter(*zip(*fitted_actual))
plt.xlabel('Predicted value')
plt.ylabel('Actual value')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()

bar_color = (0.2, 0.4, 0.6, 0.6)

# Create bar graph of the main effect A
h = 0
title = "Effect of " + titles[h]
data_m = EFFECTS[h]
data_df = np.array(dof[7])*2
data_sd = standardDeviations[h]
name_low = names[h][0]*1000
name_high = names[h][1]*1000
fig, ax = plt.subplots()
width = 1
ind = np.arange(len(data_m))
ax.bar(ind, data_m, width, yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color, align='center')
ax.set_xticks(ind+width/2)
ax.set_xticklabels([name_low, name_high], minor=False)
plt.title(title)
plt.xlabel("level of " + titles[h])
plt.ylabel('Effect with confidence interval')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()

# Create bar graph of the main effect B
h = 1
title = "Effect of " + titles[h]
data_m = EFFECTS[h]
data_df = np.array(dof[7])*2
data_sd = standardDeviations[h]
name_low = names[h][0]
name_high = names[h][1]
fig, ax = plt.subplots()
width = 1
ind = np.arange(len(data_m))
ax.bar(ind, data_m, width, yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color, align='center')
ax.set_xticks(ind+width/2)
ax.set_xticklabels([name_low, name_high], minor=False)
plt.title(title)
plt.xlabel("level of " + titles[h])
plt.ylabel('Effect with confidence interval')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()

# Create bar graph of the main effect C
h = 2
title = "Effect of " + titles[h]
data_m = EFFECTS[h]
data_df = np.array(dof[7])*2
data_sd = standardDeviations[h]
name_low = names[h][0]
name_high = names[h][1]
fig, ax = plt.subplots()
width = 1
ind = np.arange(len(data_m))
ax.bar(ind, data_m, width, yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color, align='center')
ax.set_xticks(ind+width/2)
ax.set_xticklabels([name_low, name_high], minor=False)
plt.title(title)
plt.xlabel("level of " + titles[h])
plt.ylabel('Effect with confidence interval')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()

# Create bar graph of the interaction AB
g = 0
h = 1
title = "Effect of interaction between \n" + titles[g] + " and " + titles[h]
x = []
y = []
for i in range(0, a):
    for j in range(0, b):
        bar_title = str(names[g][i]) + " x " + str(names[h][j])
        value = INTERACTIONS_AB[i][j]
        x.append(bar_title)
        y.append(value)

data_df = np.array(dof[7])*4
data_sd = [STD]*4

fig, ax = plt.subplots()
width = 1 # the width of the bars
ind = np.arange(len(y))  # the x locations for the groups
ax.bar(ind, y, width,  yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color, align='center')
ax.set_xticks(ind+width/2)
ax.set_xticklabels(x, minor=False)
plt.title(title)
plt.xlabel("level of " + titles[g] + " x level of " + titles[h])
plt.ylabel('Effect of interaction with confidence interval')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()

# Create bar graph of the interaction AC
g = 0
h = 2
title = "Effect of interaction between \n" + titles[g] + " and " + titles[h]
x = []
y = []
for i in range(0, a):
    for k in range(0, c):
        bar_title = str(names[g][i]) + " x " + str(names[h][k])
        value = INTERACTIONS_AC[i][k]
        x.append(bar_title)
        y.append(value)

data_df = np.array(dof[7])*4
data_sd = [STD]*4

fig, ax = plt.subplots()
width = 1 # the width of the bars
ind = np.arange(len(y))  # the x locations for the groups
ax.bar(ind, y, width,  yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color, align='center')
ax.set_xticks(ind+width/2)
ax.set_xticklabels(x, minor=False)
plt.title(title)
plt.xlabel("level of " + titles[g] + " x level of " + titles[h])
plt.ylabel('Effect of interaction with confidence interval')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()

# Create bar graph of the interaction BC
g = 1
h = 2
title="Effect of interaction between \n" + titles[g] + " and " + titles[h]
x = []
y = []
for j in range(0, b):
    for k in range(0, c):
        bar_title = str(names[g][j]) + " x " + str(names[h][k])
        value = INTERACTIONS_BC[j][k]
        x.append(bar_title)
        y.append(value)

data_df = np.array(dof[7])*4
data_sd = [STD]*4

fig, ax = plt.subplots()
width = 1 # the width of the bars
ind = np.arange(len(y))  # the x locations for the groups
ax.bar(ind, y, width,  yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color, align='center')
ax.set_xticks(ind+width/2)
ax.set_xticklabels(x, minor=False)
plt.title(title)
plt.xlabel("level of " + titles[g] + " x level of " + titles[h])
plt.ylabel('Effect of interaction with confidence interval')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()

# Create bar graph of the interaction ABC
title="Effect of interaction between " + titles[0] + ",\n " + titles[1] + " and\n " + titles[2]
x = []
y = []
for i in range(0, a):
    for j in range(0, b):
        for k in range(0, c):
            bar_title = str(names[0][i]) + " x \n" + str(names[1][j]) + " x\n " + str(names[2][k])
            value = INTERACTIONS_ABC[i][j][k]
            x.append(bar_title)
            y.append(value)

data_df = np.array(dof[7])*8
data_sd = [STD]*8

fig, ax = plt.subplots()
width = 1 # the width of the bars
ind = np.arange(len(y))  # the x locations for the groups
ax.bar(ind, y, width,  yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color, align='center')
ax.set_xticks(ind+width/2)
ax.set_xticklabels(x, minor=False)
plt.title(title)
plt.xlabel("level of " + titles[0] + "\n x level of " + titles[1] + "\n x level of " + titles[2])
plt.ylabel('Effect of interaction with confidence interval')
if create_output:
    plt.savefig("generated/"+title.replace("\n", "").replace(" ", "_")+".png")
    plt.show()
