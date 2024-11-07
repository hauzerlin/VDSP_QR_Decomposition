//############################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   (C) Copyright Laboratory System Integration and Silicon Implementation
//   All Right Reserved
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   2023 DSP VLSI Architecture Design
//   Fianl Project : QRD (QR Decomposition)
//   Author        : Wayne-WU (wayne.123452000@gmail.com)
//
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   File Name   : TESTBED.v
//   Module Name : TESTBED
//   Release version : v1.0
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//############################################################################
`timescale 1ns/10ps
`include "PATTERN.v"

`define INPUT_WIDTH  17		// [16:0]	in  S2.14
`define OUTPUT_WIDTH 17		// [16:0] 	out S2.14

module TESTBED();

wire clk;
wire rst_n;
wire signed	[`INPUT_WIDTH-1:0] in0, in1, in2, in3;

wire signed	[`OUTPUT_WIDTH-1:0]	out0, out1, out2, out3;

QRD U_QRD(
    .clk(clk),
    .rst_n(rst_n),
	.in0(in0),
	.in1(in1),
	.in2(in2),
	.in3(in3),
	.out0(out0),
	.out1(out1), 
	.out2(out2), 
	.out3(out3)
);

PATTERN U_PATTERN(
    .clk(clk),
    .rst_n(rst_n),
	.in0(in0),
	.in1(in1),
	.in2(in2),
	.in3(in3),
	.out0(out0),
	.out1(out1), 
	.out2(out2), 
	.out3(out3)
);

initial begin
	`ifdef RTL
		$fsdbDumpfile("QRD.fsdb");
		$fsdbDumpvars(0,"+mda");
		$fsdbDumpvars();
	`endif
	`ifdef GATE
		$sdf_annotate("QRD_SYN.sdf", U_QRD);
		$fsdbDumpfile("QRD_SYN.fsdb");
		$fsdbDumpvars();
	`endif
end

endmodule
