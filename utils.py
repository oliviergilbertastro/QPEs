import numpy as np
from math import floor, ceil

def print_table(a, header=None, title=None, space_between_columns=1, space_between_rows=0, sides=False):
    """
    Nicely print out a table

    a: array to print
    header: either an array of column names or a boolean. True -> first row will be header
    title: string, title that will be centered above the table
    space_between_columns: int, self-explanatory
    space_between_rows: int, self-explanatory
    """
    a = np.array(a, dtype=str)
    if type(header) == str or (type(header) == bool and header):
        pass
    elif (type(header) == bool and header == False) or header == None:
        header = None
    elif header != None:
        a = np.vstack((header, a))

    #Replace (None) elements with "-":
    a[a == "None"] = "-"

    #Get longest string in each column:
    column_maxes = []
    vfunc = np.vectorize(lambda x: len(x))
    a_lens = vfunc(a)
    for i in range(a.shape[1]):
        column_maxes.append(np.max(a_lens[:,i]))
    total_length = np.sum(column_maxes)+(len(column_maxes)-1)*space_between_columns #To include spaces between each column

    #Actually start printing table:
    top_and_bottom_bounds = ("─"*total_length, "─"*total_length)
    if sides:
        top_and_bottom_bounds = ("┌"+"─"*(total_length+2)+"┐", "└"+"─"*(total_length+2)+"┘")
    print()
    usable_length = total_length+4 if sides else total_length
    if title != None:
        title = floor((usable_length-len(title))/2)*" "+title+ceil((usable_length-len(title))/2)*" "
        print(f"\x1b[33m{title}\x1b[0m")
    print(top_and_bottom_bounds[0])
    #Print each row:
    for row in range(a.shape[0]):
        row_string = ""
        for column in range(a.shape[1]):
            row_string += a[row, column] + " "*(column_maxes[column]-a_lens[row,column])
            if column < a.shape[1]-1:
                row_string += " "*space_between_columns
        if row == 0 and header != None:
            row_string = f"\x1b[33m{row_string}\x1b[0m"
        if sides:
            row_string = f"│ {row_string} │"
        if row != (a.shape[0]-1):
            if sides:
                row_string += f"\n│ {' '*(total_length)} │"*space_between_rows
                if row == 0 and header != None:
                    row_string += f"\n│ {' '*(total_length)} │"
            else:
                row_string += "\n"*space_between_rows
                if row == 0 and header != None:
                    row_string += "\n"
        print(row_string)
    print(top_and_bottom_bounds[1])
    print()

if __name__ == "__main__":
    print_table([["potato", 5178, 13095], ["123", None, 1023]],
                space_between_columns=4,
                space_between_rows=0,
                header=False,
                title="My Table",
                sides=True)