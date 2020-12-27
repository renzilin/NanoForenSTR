def func_read_in_pattern_file(pattern_file_path):
    pattern_dict = {}
    with open(pattern_file_path, 'r') as file:
        for line in file:
            line_lst = line.strip().split(",")
            pattern_dict[line_lst[0]] = [line_lst[1], line_lst[2], line_lst[3], line_lst[4]]
    return pattern_dict