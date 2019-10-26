import numpy as np


def create_table(array, row_lables, column_labels, title, label):
    output = ""
    output += "\\begin{table}[H]\n"
    output += "\\caption{" + title + "} \\label{table:" + label + "}\n"
    output += "\\centering\n"
    output += "\\resizebox{\\textwidth}{!}{\\begin{tabular}{|"
    for i in range(0, len(array[0]) + 1):
        output += "l|"
    output += "}\n"
    output += "\\hline\n"

    for col in range(0, len(array[0])):
        output += column_labels[col] + " & "
    output += column_labels[len(column_labels) - 1]
    output += "\\\\ \\hline \n"

    for row in range(0, len(array)):
        output += str(row_lables[row])
        for col in range(0, len(array[0])):
                val = array[row][col]
                if isinstance(val, list) or isinstance(val, np.ndarray):
                    res = "({} , {})".format(val[0], val[1])
                elif val == 0.0:
                    res = ""
                else:
                    res = str(val)
                output += " & " + res
        output += "\\\\ \\hline \n"
    output += "\\end{tabular}}\n"
    output += "\\end{table}\n"

    print output