//############################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   (C) Copyright Laboratory System Integration and Silicon Implementation
//   All Right Reserved
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   File Name   : PATTERN.v
//   Module Name : PATTERN
//   Release version : v1.0
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//############################################################################
`ifdef RTL
	`timescale 1ns/10ps
	`include "QRD.v"
	`define CYCLE_TIME 24
`endif
`ifdef GATE
	`timescale 1ns/10ps
	`include "QRD_SYN.v"
	`define CYCLE_TIME 24
`endif
//`include "success.sv"
//`include "fail.sv"

`define INPUT_WIDTH  17		// [16:0]	in  S2.14
`define OUTPUT_WIDTH 17		// [16:0] 	out S2.14
`define PAT_NUM 8000

module PATTERN(
    clk,
    rst_n,
	in0,
	in1,
	in2,
	in3,
	out0,
	out1,
	out2,
	out3
);

output reg 	clk, rst_n;
output reg signed	[`INPUT_WIDTH-1  : 0] in0, in1, in2, in3;
input 	   signed	[`OUTPUT_WIDTH-1 : 0] out0, out1, out2, out3;
//================================================================
// wires & registers
//================================================================
reg signed [`INPUT_WIDTH-1 : 0]		input_pattern_0 	[0 : `PAT_NUM-1];
reg signed [`INPUT_WIDTH-1 : 0]		input_pattern_1 	[0 : `PAT_NUM-1];
reg signed [`INPUT_WIDTH-1 : 0]		input_pattern_2 	[0 : `PAT_NUM-1];
reg signed [`INPUT_WIDTH-1 : 0]		input_pattern_3 	[0 : `PAT_NUM-1];

reg signed [`OUTPUT_WIDTH-1 : 0]	golden_output_0 	[0 : `PAT_NUM-1];
reg signed [`OUTPUT_WIDTH-1 : 0]	golden_output_1 	[0 : `PAT_NUM-1];
reg signed [`OUTPUT_WIDTH-1 : 0]	golden_output_2 	[0 : `PAT_NUM-1];
reg signed [`OUTPUT_WIDTH-1 : 0]	golden_output_3 	[0 : `PAT_NUM-1];

reg [13:0] idx, idy;
//================================================================
// parameters & integer
//================================================================

//================================================================
// clock
//================================================================
initial	clk = 0;
always	#(`CYCLE_TIME/2.0) clk = ~clk;
//================================================================
// initial
//================================================================
// Read pattern
initial begin
	$readmemb("../00_TESTBED/in_0.txt"		, input_pattern_0);
	$readmemb("../00_TESTBED/in_1.txt"		, input_pattern_1);
	$readmemb("../00_TESTBED/in_2.txt"		, input_pattern_2);
	$readmemb("../00_TESTBED/in_3.txt"		, input_pattern_3);

	$readmemb("../00_TESTBED/out_0.txt"		 , golden_output_0);
	$readmemb("../00_TESTBED/out_1.txt"		 , golden_output_1);
	$readmemb("../00_TESTBED/out_2.txt"		 , golden_output_2);
	$readmemb("../00_TESTBED/out_3.txt"		 , golden_output_3);
end

// Simulation Start
initial	begin
	$display("-----------------------------------------------------\n");
 	$display("START!!! Simulation Start .....\n");
 	$display("-----------------------------------------------------\n");	
end

// rst_n
initial begin
		#0 rst_n = 1'b1;
		#(`CYCLE_TIME*2)  	rst_n = 1'b0;
		#(`CYCLE_TIME)  	rst_n = 1'b1;
end

// in0, in1, in2, in3;
initial begin: read_input
	integer i;
	in0 	= $signed(`OUTPUT_WIDTH'b0);
	in1 	= $signed(`OUTPUT_WIDTH'b0);
	in2 	= $signed(`OUTPUT_WIDTH'b0);
	in3 	= $signed(`OUTPUT_WIDTH'b0);
	@(negedge rst_n);
	@(posedge clk);
	
	idy = 0;
	for(i = 0; i <= `PAT_NUM-1; i = i+1)	begin	
		@(negedge clk);
		idy = idy+1;
		in0 	= input_pattern_0[i];
		in1 	= input_pattern_1[i];
		in2 	= input_pattern_2[i];
		in3 	= input_pattern_3[i];
	end
	
	@(negedge clk)		
	in0 	= $signed(`OUTPUT_WIDTH'b0);
	in1 	= $signed(`OUTPUT_WIDTH'b0);
	in2 	= $signed(`OUTPUT_WIDTH'b0);
	in3 	= $signed(`OUTPUT_WIDTH'b0);
	
	#(`CYCLE_TIME*20); $finish;
end

//================================================================
// check answer
//================================================================

initial begin: check_ans
	integer j;
	idx = 0;
	@(negedge rst_n);
	@(posedge clk);
	#(`CYCLE_TIME*6);

	for(j = 0; j <= `PAT_NUM-1; j = j+1)	begin
		//@(posedge clk);
		//idx = idx+1;
		@(negedge clk)
		if(out0!==golden_output_0[j])	begin
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		$display ("                                                                        FAIL!                                                               \n");
		$display ("                                                                Pattern NO.%03d		                                                      \n", j+1);
		$display ("                                                       Your 1's row output -> out: %d                                                		\n", out0);
		$display ("                                                     Golden 1's row output -> out: %d                                                 	\n",golden_output_0[j]);
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		end
		
		else if(out1!==golden_output_1[j])	begin
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		$display ("                                                                        FAIL!                                                               \n");
		$display ("                                                                Pattern NO.%03d		                                                      \n", j+1);
		$display ("                                                      Your 2's row output -> out: %d                                     				 	\n", out1);
		$display ("                                                    Golden 2's row output -> out: %d                                      			 \n",golden_output_1[j]);
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		end
		
		else if(out2!==golden_output_2[j])	begin
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		$display ("                                                                        FAIL!                                                               \n");
		$display ("                                                                Pattern NO.%03d		                                                      \n", j+1);
		$display ("                                                       Your 3's row output -> out: %d                                          				 \n", out2);
		$display ("                                                     Golden 3's row output -> out: %d                                           			 \n",golden_output_2[j]);
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		end

		else if(out3!==golden_output_3[j])	begin
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		$display ("                                                                        FAIL!                                                               \n");
		$display ("                                                                Pattern NO.%03d		                                                      \n", j+1);
		$display ("                                                       Your 4's row output -> out: %d                                         				 \n", out3);
		$display ("                                                     Golden 4's row output -> out: %d                                          			 \n",golden_output_3[j]);
		$display ("--------------------------------------------------------------------------------------------------------------------------------------------\n");
		end
		
		else		begin
		$display("PASS PATTERN NO.%4d", j+1);
		end
	end
end

endmodule
