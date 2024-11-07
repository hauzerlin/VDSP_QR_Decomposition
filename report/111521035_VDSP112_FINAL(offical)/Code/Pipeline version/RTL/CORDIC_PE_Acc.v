`ifndef INPUT_LENGTH
	`define INPUT_LENGTH  17
`endif

`ifndef ANGLE_LENGTH
	`define ANGLE_LENGTH  16
`endif

module CORDIC_PE_Acc(x, y, ang, mode, mu_in, map_in, rote_zero, x_out, y_out, ang_out, mu_out);

//parameter INPUT_LENGTH = 14;
//parameter ANGLE_LENGTH = 15;
parameter SHIFT_STAGE = 0;
parameter signed ELEMENTARY_ANGLE = 16'b0011001001000011;

input signed [`INPUT_LENGTH-1:0] x, y;
input signed [`ANGLE_LENGTH-1:0] ang;
input mode, mu_in, map_in, rote_zero;
output signed [`INPUT_LENGTH-1:0] x_out, y_out;
output signed [`ANGLE_LENGTH-1:0] ang_out;
output mu_out;

// caculation
wire mu;
assign mu = (rote_zero==1'b0) ? (mode==1'b0)? ((y<0)? 1'b1 : 1'b0) : ((map_in==1'b0)? mu_in : ~mu_in): ((ang<0)? 1'b0:1'b1);
assign mu_out = (rote_zero==1'b0) ? (mode==1'b0)? mu : mu_in : mu;

wire signed [`INPUT_LENGTH-1:0] Y_SHIFT, X_SHIFT;
assign X_SHIFT = (x>>>SHIFT_STAGE);
assign Y_SHIFT = (y>>>SHIFT_STAGE);

wire signed [`INPUT_LENGTH-1:0] Y_pre_add, X_pre_add;
assign X_pre_add = (mu==1'b1)? X_SHIFT : -(X_SHIFT);
assign Y_pre_add = (mu==1'b0)? Y_SHIFT : -(Y_SHIFT);
wire signed [`ANGLE_LENGTH-1:0] ang_pre_add;
assign ang_pre_add = (mu==1'b0) ? ELEMENTARY_ANGLE : -(ELEMENTARY_ANGLE);

assign x_out = x + Y_pre_add;
assign y_out = X_pre_add + y;
assign ang_out = (rote_zero==1'b0)? (x|y)? ang + ang_pre_add : $signed(0) : ang+ang_pre_add;

endmodule
