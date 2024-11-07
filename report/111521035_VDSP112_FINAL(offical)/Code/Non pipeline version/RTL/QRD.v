//############################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   2023 DSP VLSI Architecture Design
//   Fianl Project   	: QRD (QR Decomposition by HI version)
//   Author         	: Wayne-WU (wayne.123452000@gmail.com)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   File Name   : QRD.v
//   Module Name : QRD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//############################################################################

// Define data width
`define INPUT_WIDTH  17		// [16:0]	in  S2.14
`define OUTPUT_WIDTH 17		// [16:0] 	out S2.14

module QRD(
    clk,
    rst_n,
	in0, in1, in2, in3,
	out0, out1, out2, out3
);
//================================================================
//  INPUT AND OUTPUT DECLARATION                         
//================================================================
input clk;
input rst_n;
input signed	[`INPUT_WIDTH-1:0]	in0, in1, in2, in3;

output reg signed	[`OUTPUT_WIDTH-1:0]	out0, out1, out2, out3;
//================================================================
//  integer / genvar / parameters
//================================================================
parameter STAGE = 13;
parameter DIMENTION = 4;
parameter PE = 6;
integer rst_dff;
integer shift_dff;
//================================================================
//   Wires & Registers 
//================================================================
reg signed [`INPUT_WIDTH - 1:0] inreg1;				// skew input
reg signed [`INPUT_WIDTH - 1:0] inreg2 [0:1];
reg signed [`INPUT_WIDTH - 1:0] inreg3 [0:2];

reg signed [`OUTPUT_WIDTH - 1:0] delay_unit [0:DIMENTION-1];

wire signed [`OUTPUT_WIDTH - 1:0] row0 [0:DIMENTION-2];
wire signed [`OUTPUT_WIDTH - 1:0] row1 [0:DIMENTION-3];
wire signed [`OUTPUT_WIDTH - 1:0] row2;

reg signed [`OUTPUT_WIDTH - 1:0] row0_dff [0:DIMENTION-3];
reg signed [`OUTPUT_WIDTH - 1:0] row1_dff;

wire signed [`OUTPUT_WIDTH - 1:0] col0 [0:DIMENTION-2];
wire signed [`OUTPUT_WIDTH - 1:0] col1 [0:DIMENTION-3];
wire signed [`OUTPUT_WIDTH - 1:0] col2;

reg signed [`OUTPUT_WIDTH - 1:0] col0_dff [0:DIMENTION-2];
reg signed [`OUTPUT_WIDTH - 1:0] col1_dff [0:DIMENTION-3];

reg  mapping_in	[0:PE-1];
wire mapping_out[0:PE-1];
reg  signed [`OUTPUT_WIDTH - 1:0] theta_in 	[0:PE-1];
wire signed [`OUTPUT_WIDTH - 1:0] theta_out	[0:PE-1];

reg [2:0] cnt8;
reg mode [0:PE];			// mode
reg rot_zero_flag [0:2];	// Rotate 0 degree flag

// Eliminate output skew
reg signed [`OUTPUT_WIDTH - 1:0] outreg0 [0:1];
reg signed [`OUTPUT_WIDTH - 1:0] outreg1;

// Register in
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)	begin
		inreg1 		<= $signed(`INPUT_WIDTH'd0);
		inreg2[0] 	<= $signed(`INPUT_WIDTH'd0);
		inreg2[1] 	<= $signed(`INPUT_WIDTH'd0);
		inreg3[0] 	<= $signed(`INPUT_WIDTH'd0);
		inreg3[1]	<= $signed(`INPUT_WIDTH'd0);
		inreg3[2]   <= $signed(`INPUT_WIDTH'd0);
	end
	else	begin
		inreg1 		<= in1;
		inreg2[0] 	<= in2;
		inreg2[1] 	<= inreg2[0];
		inreg3[0] 	<= in3;
		inreg3[1]	<= inreg3[0];
		inreg3[2]   <= inreg3[1];
	end
end

// Delay unit
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)	begin
		for(rst_dff = 0; rst_dff <= DIMENTION-1; rst_dff = rst_dff + 1)		begin
			delay_unit[rst_dff] <= $signed(`OUTPUT_WIDTH'd0);
		end
	end
	else	begin
		delay_unit[0] <= in0;
		delay_unit[1] <= col0_dff[0];
		delay_unit[2] <= col1_dff[0];
		delay_unit[3] <= col2;
	end
end

// Row and Col register
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)	begin
		row0_dff[0] <= $signed(`OUTPUT_WIDTH'd0);
		row0_dff[1] <= $signed(`OUTPUT_WIDTH'd0);
		row1_dff	<= $signed(`OUTPUT_WIDTH'd0);
		col0_dff[0] <= $signed(`OUTPUT_WIDTH'd0);
		col0_dff[1] <= $signed(`OUTPUT_WIDTH'd0);
		col0_dff[2] <= $signed(`OUTPUT_WIDTH'd0);
		col1_dff[0]	<= $signed(`OUTPUT_WIDTH'd0);
		col1_dff[1]	<= $signed(`OUTPUT_WIDTH'd0);
	end
	else	begin
		row0_dff[0] <= row0[0];
		row0_dff[1] <= row0[1];
		row1_dff	<= row1[0];
		
		col0_dff[0] <= col0[0];
		col0_dff[1] <= col0[1];
		col0_dff[2] <= col0[2];
		col1_dff[0]	<= col1[0];
		col1_dff[1]	<= col1[1];
	end
end

// counter-8
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)				cnt8 <= 3'd0;
	else if(cnt8 == 3'd7)	cnt8 <= 3'd0;
	else					cnt8 <= cnt8 + 3'd1;
end

// mode control signal
// mode0
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)				mode[0] <= 1'd0;
	else if(cnt8 == 3'd0)	mode[0] <= 1'd1;
	else					mode[0] <= 1'd0;
end

// mode1-5
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)	begin
		for(rst_dff = 1; rst_dff <= PE; rst_dff = rst_dff + 1)		begin
			mode[rst_dff] 	<= 1'd0;
		end
	end
	else	begin
		for(shift_dff = 1; shift_dff <= PE; shift_dff = shift_dff + 1)		begin
			mode[shift_dff]	<= mode[shift_dff - 1];
		end
	end
end

// Rotate 0 degree flag
always@ (posedge clk or negedge rst_n)	begin
	if(!rst_n)				rot_zero_flag[0] <= 1'd0;
	else if(cnt8 == 3'd2)	rot_zero_flag[0] <= 1'd1;
	else					rot_zero_flag[0] <= 1'd0;
end

always@ (posedge clk or negedge rst_n)	begin
	if(!rst_n)				rot_zero_flag[1] <= 1'd0;
	else if(cnt8 == 3'd3)	rot_zero_flag[1] <= 1'd1;
	else					rot_zero_flag[1] <= 1'd0;
end

always@ (posedge clk or negedge rst_n)	begin
	if(!rst_n)								rot_zero_flag[2] <= 1'd0;
	else if((cnt8 == 3'd4)|(cnt8 == 3'd5))	rot_zero_flag[2] <= 1'd1;
	else									rot_zero_flag[2] <= 1'd0;
end

// Restore theta and mapping
genvar idx;
generate
	for(idx = 0; idx <= PE-1; idx = idx + 1)	begin:theta_mapping
		if(idx <= 4)	begin:first_and_second_row
			always@ (posedge clk or negedge rst_n)		begin
				if(!rst_n)	begin
					theta_in[idx] 	<= $signed(`OUTPUT_WIDTH'd0);
					mapping_in[idx] <= 1'd0;
				end
				else if(mode[idx])	begin
					theta_in[idx] 	<= theta_out[idx];
					mapping_in[idx] <= mapping_out[idx];
				end
			end
		end
		else if(idx <= 5)	begin:third_row
			always@ (posedge clk or negedge rst_n)		begin
				if(!rst_n)	begin
					theta_in[idx] 	<= $signed(`OUTPUT_WIDTH'd0);
					mapping_in[idx] <= 1'd0;
				end
				else if(mode[idx+1])	begin
					theta_in[idx] 	<= theta_out[idx];
					mapping_in[idx] <= mapping_out[idx];
				end
			end
		end
	end
endgenerate

// Processing Element
RPE PE0(.clk(clk), .rst_n(rst_n), .in_mode(mode[0]), .in_x(delay_unit[0]), .in_y(inreg1),
		.in_theta(theta_in[0]), .in_mapping(mapping_in[0]), .rot_zero_flag(1'd0), .out_x(row0[0]), .out_y(col0[0]),
		.out_theta(theta_out[0]), .out_mapping(mapping_out[0]));

RPE PE1(.clk(clk), .rst_n(rst_n), .in_mode(mode[1]), .in_x(row0_dff[0]), .in_y(inreg2[1]),
		.in_theta(theta_in[1]), .in_mapping(mapping_in[1]), .rot_zero_flag(1'd0), .out_x(row0[1]), .out_y(col0[1]),
		.out_theta(theta_out[1]), .out_mapping(mapping_out[1]));
		
RPE PE2(.clk(clk), .rst_n(rst_n), .in_mode(mode[2]), .in_x(row0_dff[1]), .in_y(inreg3[2]),
		.in_theta(theta_in[2]), .in_mapping(mapping_in[2]), .rot_zero_flag(1'd0), .out_x(row0[2]), .out_y(col0[2]),
		.out_theta(theta_out[2]), .out_mapping(mapping_out[2]));

RPE PE3(.clk(clk), .rst_n(rst_n), .in_mode(mode[3]), .in_x(delay_unit[1]), .in_y(col0_dff[1]),
		.in_theta(theta_in[3]), .in_mapping(mapping_in[3]), .rot_zero_flag(rot_zero_flag[0]), .out_x(row1[0]), .out_y(col1[0]),
		.out_theta(theta_out[3]), .out_mapping(mapping_out[3]));
		
RPE PE4(.clk(clk), .rst_n(rst_n), .in_mode(mode[4]), .in_x(row1_dff), .in_y(col0_dff[2]),
		.in_theta(theta_in[4]), .in_mapping(mapping_in[4]), .rot_zero_flag(rot_zero_flag[1]), .out_x(row1[1]), .out_y(col1[1]),
		.out_theta(theta_out[4]), .out_mapping(mapping_out[4]));

RPE PE5(.clk(clk), .rst_n(rst_n), .in_mode(mode[6]), .in_x(delay_unit[2]), .in_y(col1_dff[1]),
		.in_theta(theta_in[5]), .in_mapping(mapping_in[5]), .rot_zero_flag(rot_zero_flag[2]), .out_x(row2), .out_y(col2),
		.out_theta(theta_out[5]), .out_mapping(mapping_out[5]));

// Eliminate output skew
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)	begin
		outreg0[0]	<= $signed(`INPUT_WIDTH'd0);
		outreg0[1] 	<= $signed(`INPUT_WIDTH'd0);
		outreg1 	<= $signed(`INPUT_WIDTH'd0);

	end
	else	begin
		outreg0[0] 	<= row0[2];
		outreg0[1] 	<= outreg0[0];
		outreg1 	<= row1[1];
	end
end

// Register out
always@ (posedge clk or negedge rst_n)		begin
	if(!rst_n)	begin
		out0	<= $signed(`INPUT_WIDTH'd0);
		out1	<= $signed(`INPUT_WIDTH'd0);
		out2	<= $signed(`INPUT_WIDTH'd0);
		out3	<= $signed(`INPUT_WIDTH'd0);
	end
	else	begin
		out0	<= outreg0[1];
		out1	<= outreg1;
		out2	<= row2;
		out3	<= col2;
	end
end

endmodule

module RPE(
    clk,
    rst_n,
	in_mode,
	in_x,
	in_y,
	in_theta,
	in_mapping,
	rot_zero_flag,
	out_x,
	out_y,
	out_theta,
	out_mapping
);

input clk, rst_n;
input in_mode;		// mode = 0 => Rotation mode; 1 => Vectoring mode
input signed	[`INPUT_WIDTH-1:0]	in_x, in_y, in_theta;
input 	 	 						in_mapping, rot_zero_flag;

output   signed	[`OUTPUT_WIDTH-1:0] out_x, out_y, out_theta;
output 		 	 					out_mapping;
//output reg								out_mode;

// Parameter
parameter STAGE = 13;
integer rst_dff;
integer shift_dff;
//================================================================
//   Wires & Registers 
//================================================================
wire signed [`INPUT_WIDTH - 1:0] ele_angle 			[0:STAGE-1];
wire signed [`INPUT_WIDTH - 1:0] x_map, y_map, theta_map;

wire signed [`INPUT_WIDTH - 1:0] x_dff_in 			[0:STAGE-1];
wire signed [`INPUT_WIDTH - 1:0] y_dff_in 			[0:STAGE-1];
wire signed [`INPUT_WIDTH - 1:0] theta_dff_in 		[0:STAGE-1];
wire 							 control, mapping_map, mapping;

//================================================================
//   Description of circuit behavior 
//================================================================

// Elementary angle
assign ele_angle[0] = $signed(`INPUT_WIDTH'd12867);
assign ele_angle[1] = $signed(`INPUT_WIDTH'd7596);
assign ele_angle[2] = $signed(`INPUT_WIDTH'd4013);
assign ele_angle[3] = $signed(`INPUT_WIDTH'd2037);
assign ele_angle[4] = $signed(`INPUT_WIDTH'd1022);
assign ele_angle[5] = $signed(`INPUT_WIDTH'd511);
assign ele_angle[6] = $signed(`INPUT_WIDTH'd255);
assign ele_angle[7] = $signed(`INPUT_WIDTH'd127);
assign ele_angle[8] = $signed(`INPUT_WIDTH'd63);
assign ele_angle[9] = $signed(`INPUT_WIDTH'd31);
assign ele_angle[10] = $signed(`INPUT_WIDTH'd15);
assign ele_angle[11] = $signed(`INPUT_WIDTH'd7);
assign ele_angle[12] = $signed(`INPUT_WIDTH'd3);

// Initial Stage
IS IS0(.in_x(in_x), .in_y(in_y), .in_theta(in_theta), .in_mapping(in_mapping), .mode(in_mode), .rot_zero_flag(rot_zero_flag),
	   .x_map(x_map) , .y_map(y_map) , .theta_map(theta_map), .mapping_map(mapping_map));

// Choose different mapping according to different modes
assign mapping = (in_mode) ? mapping_map : in_mapping;

// Micro Rotation
genvar idx;
generate
	MR#(0) MR0(.in_x(x_map), .in_y(y_map), .in_theta(theta_map), .ele_angle(ele_angle[0]),
			   .mode(in_mode), .out_x(x_dff_in[0]), .out_y(y_dff_in[0]), .out_theta(theta_dff_in[0]));
	for(idx = 1; idx <= STAGE-1; idx = idx + 1)		begin:MR_inst
		MR#(idx) MR1(.in_x(x_dff_in[idx-1]), .in_y(y_dff_in[idx-1]), .in_theta(theta_dff_in[idx-1]), .ele_angle(ele_angle[idx]),
					 .mode(in_mode), .out_x(x_dff_in[idx]), .out_y(y_dff_in[idx]), .out_theta(theta_dff_in[idx]));
	end
endgenerate

// Scaling Factor
SF SF0(.in_x(x_dff_in[STAGE-1]), .in_y(y_dff_in[STAGE-1]), .in_theta(theta_dff_in[STAGE-1]), .in_mapping(mapping),
	   .mode(in_mode), .out_x(out_x), .out_y(out_y), .out_theta(out_theta));

assign 	out_mapping = mapping;

endmodule


module IS(in_x, in_y, in_theta, in_mapping, mode, rot_zero_flag, x_map, y_map, theta_map, mapping_map);	// Initial Stage

input  signed		[`INPUT_WIDTH-1:0]	in_x, in_y, in_theta;
input 									mode, in_mapping, rot_zero_flag;

output reg signed	[`OUTPUT_WIDTH-1:0] x_map, y_map, theta_map;
output reg				 				mapping_map;

wire  									control;

// Use mode to select control signals
assign control = (mode)? in_x[`INPUT_WIDTH-1] : in_mapping;

always @(*)	begin
	if(mode)	begin		// Vectoring mode
		if(control)	begin
			x_map = -in_x;
			y_map = in_y;
			theta_map = $signed(`OUTPUT_WIDTH'd0);
			mapping_map = 1'd1;
		end
		else	begin
			x_map = in_x;
			y_map = in_y;
			theta_map = $signed(`OUTPUT_WIDTH'd0);
			mapping_map = 1'd0;
		end
	end
	else	begin			// Rotation mode
		if(rot_zero_flag)	begin
			x_map = in_x;
			y_map = in_y;
			theta_map = $signed(`OUTPUT_WIDTH'd0);
			mapping_map = 1'd0;
		end
		else if(control)	begin
			x_map = -in_x;
			y_map = -in_y;
			theta_map = in_theta;
			mapping_map = 1'd0;
		end
		else	begin
			x_map = in_x;
			y_map = in_y;
			theta_map = in_theta;
			mapping_map = 1'd0;
		end
	end
end
endmodule


module MR(in_x, in_y, in_theta, ele_angle, mode, out_x, out_y, out_theta); 	// Micro Rotation

input  signed	[`INPUT_WIDTH-1:0]	in_x, in_y, in_theta, ele_angle;
input 								mode;
output signed	[`OUTPUT_WIDTH-1:0] out_x, out_y, out_theta;

wire signed		[`OUTPUT_WIDTH-1:0] x_shift, y_shift, angle;
wire  								control;

parameter SHIFT = 0;

// Use mode to select control signals
assign control = (mode)? ~in_y[`INPUT_WIDTH-1] : in_theta[`INPUT_WIDTH-1] ;

assign x_shift 	= (control)? -(in_x >>> SHIFT) : (in_x >>> SHIFT);
assign y_shift 	= (control)? (in_y >>> SHIFT) : -(in_y >>> SHIFT);
assign angle	= (control)? ele_angle : -ele_angle;

assign out_x = in_x + y_shift;
assign out_y = x_shift + in_y;
assign out_theta = in_theta + angle;
endmodule


module SF(in_x, in_y, in_theta, in_mapping, mode, out_x, out_y, out_theta);		// Scaling Factor

input  signed	[`INPUT_WIDTH-1:0]	in_x, in_y, in_theta;
input  								in_mapping, mode;
output signed	[`OUTPUT_WIDTH-1:0] out_x, out_y, out_theta;

assign out_x = (((in_x >>> 1) + (in_x >>> 3)) - ((in_x >>> 6) + (in_x >>> 9))) + ((in_x >>> 14)-(in_x >>> 12));
assign out_y = (((in_y >>> 1) + (in_y >>> 3)) - ((in_y >>> 6) + (in_y >>> 9))) + ((in_y >>> 14)-(in_y >>> 12));

// mapping theta for vectoring mode
assign out_theta = (mode)? ((in_mapping) ? in_theta : -in_theta) : in_theta;
endmodule