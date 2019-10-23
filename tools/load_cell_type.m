function [channel_name, cell_type] = load_cell_type(csv_filename)

table = readtable(csv_filename);

channel_name = reshape(table{:,1}, 1, []);
cell_type = reshape(table{:,2}, 1, []);


return
%%

[channel_name, cell_type] = load_cell_type('data/20180626_PSTH_CellType.csv')

%%

