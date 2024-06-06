import numpy as np
from math import floor, ceil

def print_table(a, header=None, title=None, space_between_columns=2, space_between_rows=0, borders=1, header_color="yellow", border_color="grey"):
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
    else:
        a = np.vstack((header, a))

    #Initialize the ascii characters to create the table depending on the borders parameter:
    if borders == None or borders == 0 or borders == False or borders == "none":
        characters = [" "," "," "," "," "," "]
    elif borders == "bold" or borders == 2:
        characters = ["═","║","╔","╗","╚","╝"]
    elif borders == 1 or borders == True or borders == "normal":
        characters = ["─","│","┌","┐","└","┘"]
    else:
        if type(borders) == str and len(borders) == 1:
            characters = [*(borders*6)]
        else:
            raise ValueError(f"Border style '{borders}' does not exist, use the keyword 'none', 'normal' or 'bold'.")
    
    possible_colors = ["black","red","green","yellow","blue","magenta","cyan","white"]
    #Initialize the colors:
    #Header color
    if header_color == None or header_color == "grey":
        header_color = "0"
    elif type(header_color) == str:
        header_color = header_color.lower()
        if header_color in possible_colors:
            header_color = str(possible_colors.index(header_color)+30)
        else:
            print(f"Color '{header_color}' not implemented, defaulting to grey.\nPossible colors are: {['grey']+possible_colors}")
            header_color = "0"
    else:
        raise ValueError(f"Parameter 'header_color' needs to be a string.")
    #Borders color
    if border_color == None or border_color == "grey":
        border_color = "0"
    elif type(border_color) == str:
        border_color = border_color.lower()
        if border_color in possible_colors:
            border_color = str(possible_colors.index(border_color)+30)
        else:
            print(f"Color '{border_color}' not implemented, defaulting to grey.\nPossible colors are: {['grey']+possible_colors}")
            border_color = "0"
    else:
        raise ValueError("Parameter 'border_color' needs to be a string.")

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
    data = np.array([["potato", 5178, 13095, 3151],
            ["123", None, 1023, 51515],
            ["potato", 5178, 13012495, 51515],
            ["123", 123, 1024443, 51515],
            ["potaawddwto", 5178, 13095, 51515],
            ["123", "something", 1023, 51515],
            ["potato", 5178, 13095, 51515],
            ["123", None, 1023, 51515]])
    print_table(data,
                space_between_columns=4,
                space_between_rows=0,
                header=[f"Column{i}" for i in range(data.shape[1])],
                title="My Table",
                borders="bold",
                header_color="yellow",
                border_color="blue",
                )