from openpyxl import Workbook, load_workbook, styles
from abc import ABC, abstractmethod
from datetime import date, datetime
from tkinter import filedialog
import tkinter as tk
import pandas as pd
import numpy as np
import os


class ExportDict(ABC):

    def __init__(self):

        self.export_dict = dict()
        self.init_export_dict()

    @abstractmethod
    def init_export_dict(self, *args):
        pass

    @abstractmethod
    def append_values_to_dict(self, analysis_result, *args):
        pass


def export_to_excel(

        excel_dir_path="",
        excel_name="result.xlsx",
        add_date_time=False,
        sheet_name='data'

):

    def new_func(func):

        def new_wrapper(*args, **kwargs):

            if not os.path.isdir(excel_dir_path):

                root = tk.Tk()
                root.withdraw()

                excel_dir = filedialog.askdirectory()

            else:

                excel_dir = excel_dir_path

            if add_date_time:

                today = date.today()
                now = datetime.now()
                today_str = today.strftime("%d %b")
                now_str = now.strftime("%H.%M")

                new_excel_name = "{}_{}_{}".format(today_str, now_str, excel_name)

            else:

                new_excel_name = excel_name

            if ".xlsx" not in new_excel_name:
                new_excel_name = new_excel_name + ".xlsx"

            if os.path.isdir(excel_dir):

                output_file = os.path.join(excel_dir, new_excel_name)

                if os.path.isfile(output_file):

                    writer = pd.ExcelWriter(output_file, engine="openpyxl", mode="a", if_sheet_exists="replace")

                else:

                    writer = pd.ExcelWriter(output_file, engine="openpyxl", mode="w")

                result_dict = func(*args, **kwargs)
                df = pd.DataFrame(result_dict)
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                writer.save()

        return new_wrapper

    return new_func

def write_excel_sheet(excel_path, sheet_name, data_frame: dict, first_row=2, first_column=2, overwrite="light"):

    if not os.path.isfile(str(excel_path)):

        wb = Workbook()

    else:

        wb = load_workbook(excel_path)

    if overwrite == "hard":
        if sheet_name in wb.sheetnames:
            wb.remove_sheet(wb[sheet_name])

    elif overwrite == "none":
        if sheet_name in wb.sheetnames:
            return

    if not sheet_name in wb.sheetnames:
        wb.create_sheet(sheet_name)

    sheet = wb[sheet_name]

    col = first_column
    for key in data_frame.keys():

        row = first_row
        sub_data_frame = data_frame[key]
        n_sub_element = len(sub_data_frame["unit"])

        if n_sub_element == 0:

            sheet.merge_cells(start_row=row, start_column=col, end_row=row + 1, end_column=col)
            n_sub_element = 1

        else:

            if n_sub_element > 3:
                # Add a space between data entry if n_sub_element > 3
                col += 1

            if n_sub_element > 1:
                sheet.merge_cells(start_row=row, start_column=col, end_row=row, end_column=col + n_sub_element - 1)

            for n in range(n_sub_element):

                row = first_row + 1
                cell = sheet.cell(row, col + n, value=sub_data_frame["unit"][n])
                cell.alignment = styles.Alignment(horizontal="center", vertical="center")
                cell.font = styles.Font(italic=True, size=10)

        cell = sheet.cell(first_row, col, value=key)
        cell.alignment = styles.Alignment(horizontal="center", vertical="center")
        cell.font = styles.Font(bold=True)

        for n in range(n_sub_element):

            row = first_row + 3
            data_list = (sub_data_frame["values"])[n]

            for data in data_list:

                sheet.cell(row, col, value=data)
                row += 1

            col += 1

    wb.save(excel_path)


def export_profiles_to_excel(file_path, data_input):

    time_list = data_input["time_list"]
    t_out_list = data_input["t_out_list"]
    w_out_list = data_input["w_out_list"]
    t_profile_list = data_input["t_profile_list"]
    p_profile_list = data_input["p_profile_list"]
    profile_positions = data_input["profile_positions"]

    # Write Main Data Sheet
    main_data = {

        'Time': {"unit": ["days"], "values": [time_list]},
        'T_out': {"unit": ["°C"], "values": [t_out_list]},
        'W_out': {"unit": ["kW"], "values": [w_out_list]}

    }
    write_excel_sheet(excel_path=file_path, sheet_name='Main Results', data_frame=main_data, overwrite="hard")

    profiles_list = [t_profile_list, p_profile_list]
    profiles_names = ['Temperature Profiles', 'Pressure Profiles']
    n_days = len(time_list)

    for i in range(len(profiles_list)):
        curr_profile_name = profiles_names[i]
        curr_profile_list = np.array(profiles_list[i])

        ann_data = {

            'Time': {"unit": ["days"], "values": [time_list]},
            'Annulus': {"unit": profile_positions, "values": curr_profile_list[:, 0, :].T},

        }

        write_excel_sheet(excel_path=file_path, sheet_name=curr_profile_name, data_frame=ann_data, overwrite="hard")

        tub_data = {

            'Time': {"unit": ["days"], "values": [time_list]},
            'Tubing': {"unit": profile_positions, "values": curr_profile_list[:, 1, :].T}

        }

        write_excel_sheet(excel_path=file_path, sheet_name=curr_profile_name, data_frame=tub_data, first_row=6 + n_days)
