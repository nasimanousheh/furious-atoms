def simple_lookup(np_arr, target_item):
    for i in range(len(np_arr)):
        if np_arr[i] == target_item:
            return i
    return -1