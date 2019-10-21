import numpy as np
import math
import matplotlib.pyplot as plt

from LaTeXPrinter import createTable


def square(x): return x ** 2

# Here for printing LaTeX tables later on, please ignore for now
a_title = "Degree of Phobia"
b_title = "Type of Therapy"
c_title = "Gender"
a_names = ["Mild", "Moderate", "Severe"]
b_names = ["Desensitization", "Implosion", "Insight"]
c_names = ["male", "female"]


# Layout:
#     /
#   a/
#   /_______
#   |   b
#  c|
#   |
#
# Input definition
# a = 3
# b = 3
# c = 2
# r = 3
# DATA = np.array(
#     [
#         [
#             (10, 12, 10),
#             (12, 9, 11),
#             (13, 10, 9),
#             (16, 11, 12),
#             (14, 13, 11),
#             (17, 15, 13)
#         ],
#         [
#             (15, 12, 6),
#             (12, 10, 7),
#             (14, 11, 5),
#             (17, 14, 10),
#             (18, 13, 9),
#             (16, 12, 11)
#         ],
#         [
#             (13, 11, 10),
#             (9, 7, 6),
#             (11, 8, 8),
#             (16, 10, 11),
#             (12, 12, 10),
#             (14, 14, 9)
#         ]
#     ]
# )

#Sample data for 2^k analysis, taken from the slides of lecture 1_2 slide 31.
a = 2
b = 2
c = 2
r = 3 # number of repetitions
DATA = np.array(
    [
        [
            (86, 58),
            (80, 62),
            (74, 60),
            (34, 22),
            (30, 18),
            (35, 20)
        ],
        [
            (50, 46),
            (55, 42),
            (54, 44),
            (11, 14),
            (15, 16),
            (19, 12)
        ]
    ]
)

# # Calculate averages
AVERAGES = np.zeros((a, c, b))

for i in range(0, a):
    for k in range(0, b):
        for j in range(0, c):
            average = np.average(DATA[i][j*r : (j+1)*r], axis=0)
            AVERAGES[i][j] = average

print AVERAGES

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

print SUMS
print MEANS

# # Calculate total sum and mean
TOTAL_SUM = np.sum(AVERAGES)
TOTAL_MEAN = np.average(AVERAGES)

EFFECTS = MEANS - TOTAL_MEAN

print TOTAL_MEAN
print EFFECTS

SSA = c*b*r*np.sum(square(EFFECTS[0]))
SSB = c*a*r*np.sum(square(EFFECTS[1]))
SSC = a*b*r*np.sum(square(EFFECTS[2]))
print SSA
print SSB
print SSC

INTERACTIONS = np.zeros((a, c, b))

for C in range(0, c):
    for B in range(0, b):
        for A in range(0, a):
            INTERACTIONS[A][B][C] = AVERAGES[A][B][C] - MEANS[0][A] - MEANS[1][B] - MEANS[2][C] + TOTAL_MEAN

print INTERACTIONS

## This is the part where I'm stuck right now
SSAB = INTERACTIONS[0][0]**2 + EFFECTS[0][1]**2


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






