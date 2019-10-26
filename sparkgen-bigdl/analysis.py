import numpy as np
from scipy.stats import f, t
import math
from LaTeXPrinter import create_table
from dataCollector import collector
import matplotlib.pyplot as plt

create_output = True

def square(x): return x ** 2


titles, names, parameters, DATA = collector()

# DATA = np.log10(DATA)
# print DATA

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

# print SUMS
# print MEANS

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

# print "------------------------------------------------------------------"
SSY = np.sum(np.sum(np.sum(square(DATA))))
SSO = a * b * c * r * (TOTAL_MEAN ** 2)
SST = SSY - SSO
SS[7] = SST - np.sum(SS)

# print SSY
# print SSO
# print SST

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

# print SS

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

# print MS

# Calculate F-comp values
Fcomp = np.zeros((7, 1))
Fcomp[0] = MS[0]/MS[7]
Fcomp[1] = MS[1]/MS[7]
Fcomp[2] = MS[2]/MS[7]
Fcomp[3] = MS[3]/MS[7]
Fcomp[4] = MS[4]/MS[7]
Fcomp[5] = MS[5]/MS[7]
Fcomp[6] = MS[6]/MS[7]

# print Fcomp

alpha = 0.10
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


# print standardDeviations

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

# print CIs_MAIN

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

average_table_data = np.zeros((a * c + a, b))
average_table_data[1, :] = AVERAGES[0, 0, :]
average_table_data[2, :] = AVERAGES[0, 1, :]
average_table_data[4, :] = AVERAGES[1, 0, :]
average_table_data[5, :] = AVERAGES[1, 1, :]

# print average_table_data

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

effects_table_data = np.zeros((a + b + c + 3, 2))
effects_table_data[1:3, 0] = EFFECTS[0]
effects_table_data[1:3, 1] = standardDeviations[0]
effects_table_data[4:6, 0] = EFFECTS[1]
effects_table_data[4:6, 1] = standardDeviations[1]
effects_table_data[7:10, 0] = EFFECTS[2]
effects_table_data[7:10, 1] = standardDeviations[2]

STD_mu = s_e*math.sqrt(1.0 / (a*b*r*c))
CI_mu = [TOTAL_MEAN + t.ppf(1 - alpha, dof[7])*STD_mu, TOTAL_MEAN - t.ppf(1 - alpha, dof[7])*STD_mu]

# Create LaTex table for the confidence intervals of the effects.
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
    plt.savefig("generated/"+title+".png")
    plt.show()

bar_color = (0.2, 0.4, 0.6, 0.6)

h = 0
title = "Confidence intervals for " + titles[h]
data_m = EFFECTS[h]
data_df = np.array(dof[7])*2
data_sd = standardDeviations[h]
name_low = names[h][0]*1000
name_high = names[h][1]*1000
plt.bar([name_low,name_high], data_m, yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color)
plt.title(title)
if create_output:
    plt.savefig("generated/"+title+".png")
    plt.show()

h = 1
title = "Confidence intervals for " + titles[h]
data_m = EFFECTS[h]
data_df = np.array(dof[7])*2
data_sd = standardDeviations[h]
name_low = names[h][0]
name_high = names[h][1]
plt.bar([name_low,name_high], data_m, yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color)
plt.title(title)
if create_output:
    plt.savefig("generated/"+title+".png")
    plt.show()

h = 2
title = "Confidence intervals for " + titles[h]
data_m = EFFECTS[h]
data_df = np.array(dof[7])*2
data_sd = standardDeviations[h]
name_low = names[h][0]
name_high = names[h][1]
plt.bar([name_low,name_high], data_m, yerr=t.ppf(1-alpha, data_df)*data_sd, color=bar_color)
plt.title(title)
if create_output:
    plt.savefig("generated/"+title+".png")
    plt.show()

g = 0
h = 1
title = "Effect of interaction between " + titles[g] + " and " + titles[h]
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
    plt.savefig("generated/"+title+".png")
    plt.show()

g = 0
h = 2
title = "Effect of interaction between " + titles[g] + " and " + titles[h]
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
    plt.savefig("generated/"+title+".png")
    plt.show()

g = 1
h = 2
title="Effect of interaction between " + titles[g] + " and " + titles[h]
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
    plt.savefig("generated/"+title+".png")
    plt.show()

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
    plt.savefig("generated/"+title.replace("\n","")+".png")
    plt.show()


# # Calculate row and column sums and means.
# ROW_SUMS = np.sum(AVERAGES, 1)
# COL_SUMS = np.sum(AVERAGES, 0)
#
# ROW_MEANS = ROW_SUMS / a
# COL_MEANS = COL_SUMS / b
#
# # Calculate total sum and mean
# TOTAL_SUM = sum(ROW_SUMS)
# TOTAL_MEAN = TOTAL_SUM / (b * a)
#
#
# # Calculate row and column effects
# ROW_EFFECTS = ROW_MEANS - TOTAL_MEAN
# COL_EFFECTS = COL_MEANS - TOTAL_MEAN
#
# # Calculate interactions
# INTERACTIONS = np.zeros((b, a))
#
# for dataset in range(0, b):
#     for algorithm in range(0, a):
#         INTERACTIONS[dataset][algorithm] = AVERAGES[dataset][algorithm] - ROW_MEANS[dataset] - COL_MEANS[
#             algorithm] + TOTAL_MEAN
#
# # Create LaTex table of the interactions
# title = "Effects of the interactions"
# label = "effects_interactions"
# createTable(INTERACTIONS, b_names,  [b_title] + a_names, title, label)
#
# # Create LaTex table of the computation of effects
# print_data = np.append(AVERAGES, ROW_SUMS.reshape((b, 1)), axis=1)
# print_data = np.append(print_data, ROW_MEANS.reshape((b, 1)), axis=1)
# print_data = np.append(print_data, ROW_EFFECTS.reshape((b, 1)), axis=1)
#
# print_data = np.append(print_data, np.append(COL_SUMS, [TOTAL_SUM, 0, 0]).reshape((1, a + 3)), axis=0)
# print_data = np.append(print_data, np.append(COL_MEANS, [0, TOTAL_MEAN, 0]).reshape((1, a + 3)), axis=0)
# print_data = np.append(print_data, np.append(COL_EFFECTS, [0, 0, 0]).reshape((1, a + 3)), axis=0)
#
# title = "Computation of effects"
# label = "computation_effects"
# createTable(print_data, b_names + ["Col Sum", "Col Mean", "Col Effect"],  [b_title] + a_names + ["Row Sum", "Row Mean", "Row Effect"], title, label)
#
# """The created ANOVA table has the following structure:
# columns:
# [X][0] Sum of squares,
# [X][1] Percentage of variation explained,
# [X][2] Degrees of freedom,
# [X][3] Mean squares,
# [X][4] F-comp. values,
# [X][5] F-table values
#
# rows:
# [0][X] y
# [1][X] bar{y}...
# [2][X] y - bar{y}...
# [3][X] Classifiers (factor A)
# [4][X] Datasets (factor B)
# [5][X] Interactions
# [6][X] Error
# """
# # Start computation of ANOVA results.
# RESULTS = np.zeros((7, 6))
#
# # Calculate  sum of squares
# RESULTS[0][0] = sum(sum(square(DATA)))
# RESULTS[1][0] = a * b * r * TOTAL_MEAN ** 2
# RESULTS[2][0] = RESULTS[0][0] - RESULTS[1][0]
# RESULTS[3][0] = b * r * sum(square(COL_EFFECTS))
# RESULTS[4][0] = a * r * sum(square(ROW_EFFECTS))
# RESULTS[5][0] = r * sum(sum(square(INTERACTIONS)))
# RESULTS[6][0] = RESULTS[2][0] - RESULTS[3][0] - RESULTS[4][0] - RESULTS[5][0]
#
# # Calculate amount of variation explained
# RESULTS[2][1] = 100
# RESULTS[3][1] = 100*(RESULTS[3][0]/RESULTS[2][0])
# RESULTS[4][1] = 100*(RESULTS[4][0]/RESULTS[2][0])
# RESULTS[5][1] = 100*(RESULTS[5][0]/RESULTS[2][0])
# RESULTS[6][1] = 100*(RESULTS[6][0]/RESULTS[2][0])
#
# # Calculate degrees of freedom
# RESULTS[0][2] = a*b*r
# RESULTS[1][2] = 1
# RESULTS[2][2] = a*b*r - 1
# RESULTS[3][2] = a - 1
# RESULTS[4][2] = b - 1
# RESULTS[5][2] = (a - 1)*(b - 1)
# RESULTS[6][2] = a*b*(r - 1)
#
# # Calculate mean squares
# RESULTS[3][3] = RESULTS[3][0] / (a - 1)
# RESULTS[4][3] = RESULTS[4][0] / (b - 1)
# RESULTS[5][3] = RESULTS[5][0] / ((a - 1)*(b - 1))
# RESULTS[6][3] = RESULTS[6][0] / (a*b*(r - 1))
#
# # Calculate F-comp values
# RESULTS[3][4] = RESULTS[3][3] / RESULTS[6][3]
# RESULTS[4][4] = RESULTS[4][3] / RESULTS[6][3]
# RESULTS[5][4] = RESULTS[5][3] / RESULTS[6][3]
#
# # Placeholders for F-table values. Could not get automatic F-table lookup to work.
# RESULTS[3][5] = 2.49
# RESULTS[4][5] = 2.14
# RESULTS[5][5] = 1.88
#
# #  Create LaTeX table of the ANOVA results.
# table_row_names = ["$y$", "$\\bar{y}...$", "$y - \\bar{y}...$"] + [a_title] + [b_title] + ["Interaction"] + ["error"]
# table_col_names = ["Component", "Sum of squares", "\\% Variation", "DF", "Mean Square", "F-Comp.", "F-Table"]
#
# title = "ANOVA results"
# label = "anova"
# createTable(RESULTS, table_row_names,  table_col_names, title, label)
#
# Se = math.sqrt(RESULTS[6][3])
#
# # Calculate standard deviations
# STD_A = Se*math.sqrt((a - 1.0) / (a*b*r))
# STD_B = Se*math.sqrt((b - 1.0) / (a*b*r))
# STD_mu = Se*math.sqrt(1.0 / (a*b*r))
#
# # t[0.95, ab(r-1)], looked up manually as I couldn't get the automatic lookup working.
# quantile = 1.697
#
# # Calculate confidence intervals
# CI_A = np.zeros((a, 3))
# CI_B = np.zeros((b, 3))
#
# CI_mu = (TOTAL_MEAN - STD_mu*quantile, TOTAL_MEAN + STD_mu*quantile)
#
# for i in range(0, a):
#     CI_A[i][0] = STD_A
#
# for i in range(0, b):
#     CI_B[i][0] = STD_B
#
# # Calculate confidence intervals for factor A
# for i in range(0, a):
#     minVal = COL_EFFECTS[i] - CI_A[i][0]*quantile
#     maxVal = COL_EFFECTS[i] + CI_A[i][0]*quantile
#     CI_A[i][1] = minVal
#     CI_A[i][2] = maxVal
#
# # Calculate confidence intervals for factor B
# for i in range(0, b):
#     minVal = ROW_EFFECTS[i] - CI_B[i][0]*quantile
#     maxVal = ROW_EFFECTS[i] + CI_B[i][0]*quantile
#     CI_B[i][1] = minVal
#     CI_B[i][2] = maxVal
#
# # Create LaTex table for the confidence intervals of the effects.
# print_data = [[TOTAL_MEAN, STD_mu, "({} , {})".format(CI_mu[0], CI_mu[1])]]
#
# print_data = np.append(print_data, [["", "", ""]], axis=0)
# for i in range(0, len(CI_A)):
#     print_data = np.append(print_data, [[COL_EFFECTS[i], STD_A, "({} , {})".format(CI_A[i][1], CI_A[i][2])]], axis=0)
#
# print_data = np.append(print_data, [["", "", ""]], axis=0)
# for i in range(0, len(CI_B)):
#     print_data = np.append(print_data, [[ROW_EFFECTS[i], STD_B, "({} , {})".format(CI_B[i][1], CI_B[i][2])]], axis=0)
#
# table_row_names = ["$\\mu$"] + [a_title] + a_names + [b_title] + b_names
# table_col_names = ["Parameter", "Mean Effect", "Standard Deviation", "CI"]
#
# title = "Confidence intervals of effects"
# label = "CI_effects"
# createTable(print_data, table_row_names, table_col_names, title, label)
#
# # Calculate confidence intervals for the interactions.
# STD_AB = Se*math.sqrt(((a - 1.0)*(b - 1)) / (a*b*r))
#
# CI_AB = np.zeros((b, a), dtype=(float, 2))
#
# for i in range(0, a):
#     for j in range(0, b):
#         minVal = INTERACTIONS[j][i] - STD_AB*quantile
#         maxVal = INTERACTIONS[j][i] + STD_AB*quantile
#         CI_AB[j][i] = (minVal, maxVal)
#
# # Create LaTex table for the confidence intervals of interactions
# table_row_names = b_names
# table_col_names = [b_title] + a_names
#
# title = "Confidence intervals of interactions"
# label = "CI_interactions_effects"
# createTable(CI_AB, table_row_names, table_col_names, title, label)
#
# # Create plots of fitted value against residual and fitted value against actual value.
# fitted_residuals = []
# fitted_actual = []
# for row in range(0, b):
#     block = DATA[row*r:(row+1)*r]
#     for block_row in block:
#         for col in range(0, a):
#             actual = block_row[col]
#             fitted = TOTAL_MEAN + ROW_EFFECTS[row] + COL_EFFECTS[col] + INTERACTIONS[row][col]
#             residual = actual - fitted
#
#             fitted_residuals.append((fitted, residual))
#             fitted_actual.append((fitted, actual))
#
# plt.figure()
# plt.subplot(211)
# plt.scatter(*zip(*fitted_residuals))
# plt.xlabel('Predicted value')
# plt.ylabel('Residual')
#
# plt.subplot(212)
# plt.scatter(*zip(*fitted_actual))
# plt.xlabel('Predicted value')
# plt.ylabel('Actual value')
# plt.show()






