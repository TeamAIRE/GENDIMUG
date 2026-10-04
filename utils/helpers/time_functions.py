from time import time, strftime, gmtime


def get_current_time():
    """ Returns current time (float)"""
    return time()


def total_time(ref_time):
    final_time = strftime("%H:%M:%S", gmtime(get_current_time() - ref_time))
    print(f"\n--- Completion time >>> {final_time} ---")
