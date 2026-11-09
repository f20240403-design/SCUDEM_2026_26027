"""Print sheets, headers and sample rows of an xlsx (data only, no formulas executed)."""
import sys
import openpyxl
sys.stdout.reconfigure(encoding='utf-8')
wb = openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
nrows = int(sys.argv[2]) if len(sys.argv) > 2 else 8
for ws in wb.worksheets:
    rows = list(ws.iter_rows(values_only=True))
    print(f'=== {ws.title}: {len(rows)} rows x {max((len(r) for r in rows), default=0)} cols')
    for r in rows[:nrows]:
        print('  ', [c for c in r][:20])
