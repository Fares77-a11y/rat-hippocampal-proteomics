import pandas as pd

files = [
    'UCP1_Rat_Hippo_Proteomics.xlsx',
    'Rat_Serum_Training_Anonymized.xlsx',
    'Rat_Serum_Pre_Diabetes_HFD_Statistical.xlsx'
]

for f in files:
    print('=' * 70)
    print(f'FILE: {f}')
    print('=' * 70)
    try:
        xl = pd.ExcelFile(f)
        print(f'  Sheets: {xl.sheet_names}')
        for sheet in xl.sheet_names:
            df = pd.read_excel(f, sheet_name=sheet, nrows=5)
            print(f'\n  --- Sheet: "{sheet}" ---')
            print(f'  Shape (preview): {df.shape}')
            print(f'  Columns: {df.columns.tolist()[:15]}')
            if len(df.columns) > 15:
                print(f'  ... and {len(df.columns) - 15} more columns')
    except Exception as e:
        print(f"Error reading {f}: {e}")
    print()
