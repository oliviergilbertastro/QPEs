import numpy as np
from math import floor, ceil

def print_table(a, header=None, title=None, space_between_columns=1, space_between_rows=0, borders=1, header_color="yellow", border_color="grey"):
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

    #Initialize the ascii characters to create the table depending on the borders parameter:
    if borders == None or borders == 0 or borders == False or borders == "none":
        characters = [" "," "," "," "," "," "]
    elif borders == "bold" or borders == 2:
        characters = ["═","║","╔","╗","╚","╝"]
    elif borders == 1 or borders == True or borders == "normal":
        characters = ["─","│","┌","┐","└","┘"]
    else:
        raise ValueError(f"Border style '{borders}' does not exist, use the keyword 'none', 'normal' or 'bold'.")
    
    #Initialize the colors:
    if header_color == None or header_color == "grey":
        header_color = "0"
    elif type(header_color) == str:
        header_color = header_color.lower()
    else:
        raise ValueError("Parameter 'header_color' needs to be a string.")
    if header_color == "black":
        header_color = "30"
    elif header_color == "red":
        header_color = "31"
    elif header_color == "green":
        header_color = "32"
    elif header_color == "yellow":
        header_color = "33"
    elif header_color == "blue":
        header_color = "34"
    elif header_color == "magenta":
        header_color = "35"
    elif header_color == "cyan":
        header_color = "36"
    elif header_color == "white":
        header_color = "37"
    else:
        header_color = "0"

    if border_color == None or border_color == "grey":
        border_color = "0"
    elif type(border_color) == str:
        border_color = border_color.lower()
    else:
        raise ValueError("Parameter 'border_color' needs to be a string.")
    if border_color == "black":
        border_color = "30"
    elif border_color == "red":
        border_color = "31"
    elif border_color == "green":
        border_color = "32"
    elif border_color == "yellow":
        border_color = "33"
    elif border_color == "blue":
        border_color = "34"
    elif border_color == "magenta":
        border_color = "35"
    elif border_color == "cyan":
        border_color = "36"
    elif border_color == "white":
        border_color = "37"
    else:
        border_color = "0"
    for i in range(len(characters)):
        characters[i] = f"\x1b[{border_color}m{characters[i]}\x1b[0m"

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
    top_and_bottom_bounds = (characters[2]+characters[0]*(total_length+2)+characters[3], characters[4]+characters[0]*(total_length+2)+characters[5])
    print()
    usable_length = total_length+4
    if title != None:
        title = floor((usable_length-len(title))/2)*" "+title+ceil((usable_length-len(title))/2)*" "
        print(f"\x1b[{header_color}m{title}\x1b[0m")
    print(top_and_bottom_bounds[0])
    #Print each row:
    for row in range(a.shape[0]):
        row_string = ""
        for column in range(a.shape[1]):
            row_string += a[row, column] + " "*(column_maxes[column]-a_lens[row,column])
            if column < a.shape[1]-1:
                row_string += " "*space_between_columns
        if row == 0 and header != None:
            row_string = f"\x1b[{header_color}m{row_string}\x1b[0m"
        row_string = f"{characters[1]} {row_string} {characters[1]}"
        if row != (a.shape[0]-1):
            row_string += f"\n{characters[1]} {' '*(total_length)} {characters[1]}"*space_between_rows
            if row == 0 and header != None:
                row_string += f"\n{characters[1]} {' '*(total_length)} {characters[1]}"
        print(row_string)
    print(top_and_bottom_bounds[1])
    print()

if __name__ == "__main__":
    print_table([["potato", 5178, 13095], ["123", None, 1023],["potato", 5178, 13012495], ["123", 123, 1024443],["potaawddwto", 5178, 13095], ["123", "something", 1023],["potato", 5178, 13095], ["123", None, 1023]],
                space_between_columns=4,
                space_between_rows=0,
                header=["Column1", "Column2", "Column3"],
                title="My Table",
                borders=2,
                header_color="yellow",
                border_color="white",
                )