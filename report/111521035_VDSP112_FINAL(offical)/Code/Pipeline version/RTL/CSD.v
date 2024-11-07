`ifndef INPUT_LENGTH
`define INPUT_LENGTH 17
`endif

module CSD(in, out);

input signed [`INPUT_LENGTH-1:0] in;
output signed [`INPUT_LENGTH-1:0] out;

wire signed [`INPUT_LENGTH-1:0] xbuff1, xbuff3, xbuff6, xbuff9, xbuff12, xbuff14;
wire signed [`INPUT_LENGTH-1:0] add_buff1, add_buff2, add_buff3, add_buff4;

assign xbuff1 = in >>> 1;
assign xbuff3 = in >>> 3;
assign xbuff6 = in >>> 6;
assign xbuff9 = in >>> 9;
assign xbuff12= in >>> 12;
assign xbuff14= in >>> 14;

assign add_buff1 = xbuff1+xbuff3;
assign add_buff2 = xbuff6+xbuff9;
assign add_buff3 = xbuff12-xbuff14;
assign out = add_buff1 - add_buff2 -add_buff3 ;

endmodule
