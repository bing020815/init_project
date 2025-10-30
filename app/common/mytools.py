#%%
import os
import pandas as pd
import json
from decimal import Decimal, ROUND_HALF_UP

## check directory and change directory
# print parent directory
parent_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
print(f"parent directory: {parent_directory}")

# Get the grandparent directory (the parent of the parent)
grandparent_directory = os.path.abspath(os.path.join(parent_directory, os.pardir))
print(f"grandparent directory: {grandparent_directory}")

# Change the current working directory 
current_directory = os.getcwd()
print(f"initial directory: {current_directory}")

# os.chdir('../raw_data')  # Change to the raw_data directory
# current_directory = os.getcwd()
# print(f"reset directory: {current_directory}")


### Time related functions
def get_unix_time():
    """
        get UTC+8 unix time
    """
    from datetime import datetime
    import pytz

    # Get the current time in UTC
    now_utc = datetime.now(pytz.timezone('UTC'))

    # Convert to UTC+8
    now_utc8 = now_utc.astimezone(pytz.timezone('Asia/Taipei'))

    # Return the Unix timestamp
    utc8_timestamp = int(now_utc8.timestamp())

    return utc8_timestamp

def get_date_time():
    """
        get UTC+8 date time
    """
    from datetime import datetime
    import pytz

    # Get the current time in UTC
    now_utc = datetime.now(pytz.timezone('UTC'))

    # Convert to UTC+8
    now_utc8 = now_utc.astimezone(pytz.timezone('Asia/Taipei'))
    now_utc8_str = now_utc8.strftime("%Y-%m-%d")

    # Return the Unix timestamp
    return now_utc8_str

def unit_time_to_date_time(unix_time):
    """
        convert unix time to date time
    """
    import datetime
    import pytz
    date_time = datetime.datetime.fromtimestamp(unix_time,
                                pytz.timezone('Asia/Taipei'))
    date_time_str = date_time.strftime("%Y-%m-%d")
    return date_time_str

def decompose_date_time(date_time_str):
    """
        decompose date time str to year, month, day
    """
    year, month, day = date_time_str.split('-')
    return year, month, day

def compare_date(target_date, date_list):
    """
        Compare target_date with date_list,
        return the date in date_list which is earlier than target_date
    """
    from datetime import datetime
    early_date = []
    target_date = datetime.strptime(target_date, '%Y-%m-%d')
    for d in date_list:
        if datetime.strptime(d, '%Y-%m-%d') < target_date:
            early_date.append(d)
    return sorted(early_date, reverse=False)

def check_date_format(df, date_column):
    """
        Check if a date column in a DataFrame has the format "YYYY-MM-DD"
    """
    try:
        # Try to convert the date column to datetime
        pd.to_datetime(df[date_column], format="%Y-%m-%d")
        # print("The date column has the format 'YYYY-MM-DD'.")
    except ValueError:
        print("The date column doesn't have the format 'YYYY-MM-DD'.")

def get_random_string(length):
    """
        get random string of length range from 0-9, a-z and A-Z
    """
    import random
    import string
    letters = string.ascii_lowercase + string.ascii_uppercase + string.digits
    result_str = ''.join(random.choice(letters) for _ in range(length))
    return result_str

def convert_datetime_to_date_str(series):
    """
        Convert datetime columns to date string format (YYYY-MM-DD)
        Only processes datetime columns or object columns with 'timestamp' in name
        
        Args:
            series: pandas Series to process
            
        Returns:
            pandas Series: Original series or converted date strings
    """
    if (series.dtype.name.startswith('datetime') or 
        (series.dtype == 'object' and 'timestamp' in series.name.lower())):
        return pd.to_datetime(series).dt.strftime('%Y-%m-%d')
    else:
        return series

def convert_cell_phone_to_str(cell_phone):
    """
    Convert a cell phone number to a standardized 10-digit string format.
    This function takes a cell phone number (which may be in various formats such as 
    int, float, or string) and converts it to a standardized 10-digit string format,
    padding with leading zeros if necessary. If the input is invalid (NaN, empty, or None),
    it returns a default string of 10 zeros.
    Args:
        cell_phone: The cell phone number to convert. Can be int, float, string, or None/NaN.
    Returns:
        str: A 10-digit string representation of the cell phone number. Returns '0000000000'
             if the input is invalid, empty, or NaN.
    Example:
        >>> convert_cell_phone_to_str(1234567890)
        '1234567890'
        >>> convert_cell_phone_to_str(123456789.0)
        '0123456789'
        >>> convert_cell_phone_to_str(None)
        '0000000000'
        >>> convert_cell_phone_to_str('')
        '0000000000'
    """
    if pd.notna(cell_phone) and str(cell_phone).strip() != '':
        return str(int(float(cell_phone))).zfill(10)
    else:
        return '0000000000'

### File related functions
def get_batch_df(file_name) -> pd.DataFrame:
    """
        get the batch data from file_name
        example:
            file_name = "patient_batch.csv",
            batch_file_path = "../raw_data/patient_batch.csv"
    """
    batch_file_path = f"../raw_data/{file_name}"
    print(f"load file from {batch_file_path}")
    # Read the CSV file
    df = pd.read_csv(batch_file_path, sep=',')
    # Check the format of the date column
    check_date_format(df, "sample_date")
    # print number of records
    print(f"number of records: {len(df)}")
    return df


def get_json_data(file_path_of_json) -> dict:
    """
        get the json data from file_path_of_json
        example:
            file_path_of_json = "../raw_data/items.json",
    """
    file_path = file_path_of_json
    print(f"load file from {file_path}")
    # Read the JSON file
    with open(file_path, "r", encoding="utf-8") as file:
        json_data = json.load(file)
    return json_data


def save_json_data(file_path_of_json, json_data) -> None:
    """
        save the json data to file_path_of_json
        example:
            file_path_of_json = "../raw_data/items.json",
    """
    file_path = file_path_of_json
    print(f"save file to {file_path}")
    # Write the JSON data to a file
    with open(file_path, "w", encoding="utf-8") as file:
        json.dump(json_data, file, ensure_ascii=False, indent=4)
    print(f"save to {file_path} done!")


# %%
# tool functions
def round_ex(num, decimal):
    str_deci = 1
    for _ in range(decimal):
        str_deci = str_deci / 10
    str_deci = str(str_deci)
    result = Decimal(str(num)).quantize(Decimal(str_deci), rounding=ROUND_HALF_UP)
    result = float(result)
    return result

