`timescale 1ns/1ps
module testbench();

//parameter cycle = 3.35;   
parameter cycle = 11.0;     
//parameter cycle =   20;   

reg clk, rst_n, valid;

parameter HOLD = 4;
parameter INPUT_LENGTH = 17;
parameter ANGLE_LENGTH = 16;

//parameter INPUT_PATTERN = 10000;
//parameter INPUT_PATTERN = 24;
parameter INPUT_PATTERN = 8000;

reg signed [INPUT_LENGTH-1:0] in1, in2, in3, in4;
//reg signed [ANGLE_LENGTH-1:0] ang_in;
//reg flag_in;
//reg mode;

wire signed [INPUT_LENGTH-1:0] out1, out2, out3, out4;
//wire signed [ANGLE_LENGTH-1:0] ang_out;
//wire flag_out;
//wire mode_out;

integer idx = 0, idy=0, err=0;
integer temp =0;

reg signed [INPUT_LENGTH-1:0] in1_data [INPUT_PATTERN-1:0], in2_data [INPUT_PATTERN-1:0], in3_data [INPUT_PATTERN-1:0], in4_data [INPUT_PATTERN-1:0];
reg signed [INPUT_LENGTH-1:0] out1_data [INPUT_PATTERN-1:0], out2_data [INPUT_PATTERN-1:0], out3_data [INPUT_PATTERN-1:0], out4_data [INPUT_PATTERN-1:0];

//CORDIC test(	.clk(clk), .rst_n(rst_n), .mode(mode), .x(x_in), .y(y_in), .angle(ang_in), .map_flag(flag_in),
//		.x_out(x_out), . y_out(y_out), .angle_out(ang_out), .map_flag_out(flag_out), .mode_out(mode_out) );
QR_decomposition_v2 test (	.clk(clk), .rst_n(rst_n), .valid(valid), 
			.in_1(in1), .in_2(in2), .in_3(in3), .in_4(in4),
			.out_1(out1), .out_2(out2), .out_3(out3), .out_4(out4));

initial begin
	$readmemb("./txt/in1.txt",in1_data);
	$readmemb("./txt/in2.txt",in2_data);
	$readmemb("./txt/in3.txt",in3_data);
	$readmemb("./txt/in4.txt",in4_data);
	$readmemb("./ans/out_0.txt",out1_data);
	$readmemb("./ans/out_1.txt",out2_data);
	$readmemb("./ans/out_2.txt",out3_data);
	$readmemb("./ans/out_3.txt",out4_data);
end

initial begin
	`ifdef RTL
	$fsdbDumpfile("RTL_PE.fsdb");
	$fsdbDumpvars(0,"+mda");
	`endif
	`ifdef GATE
	$sdf_annotate("../SYN/QR_decomposition_v2_syn.sdf",test);
	$fsdbDumpfile("SYN.fsdb");
	$fsdbDumpvars(0,"+mda");
	$fsdbDumpvars();
	`endif
end

always #(cycle/2) clk = ~clk;

initial begin
	clk = 1'b0;
	rst_n = 1'b1;
#(cycle * 1) rst_n = 1'b0;
#(cycle * 2.75) rst_n = 1'b1;
end

initial begin
err =0;
@(posedge rst_n)
#(cycle * HOLD * 5)
	for(idy=0;idy<INPUT_PATTERN;idy=idy+1)begin
		@(negedge clk)
		if(out1 != out1_data[idy])begin
			temp = 1;
			err = err +1;
			$display($time," out1 error :%4d error at %d\n", err, idy);
		end
		if(out2 != out2_data[idy])begin
			temp = 1;
			err = err +1;
			$display($time," out2 error :%4d error at %d\n", err, idy);
		end
		if(out3 != out3_data[idy])begin
			temp = 1;
			err = err +1;
			$display($time," out3 error :%4d error at %d\n", err, idy);
		end
		if(out4 != out4_data[idy])begin
			temp = 1;
			err = err +1;
			$display($time," out4 error :%4d error at %d\n", err, idy);
		end
		if(temp == 0)begin
			$display("NO.%4d pattern is success!\n",idy); 
			temp =0;
		end
		temp =0;
		/*if(temp == 0)begin
			$display("NO.%4d pattern is success!\n",idy); 
			temp =0;
		end
		temp =0;*/
	end
	if(err == 0) begin
		$display("***********************************\n");
		$display(" The simulation is successful ! \n");
		$display("***********************************\n");
	end
	else begin
		$display("***********************************\n");
		$display(" The simulation is fail !\n");
		$display(" There is %4d error\n", err);
		$display("***********************************\n");
	end
end
initial begin
idx = 0;
valid = 0;
in1 = 0;
in2 = 0;
in3 = 0;
in4 = 0;
@(posedge rst_n)
	for (idx=0;idx<INPUT_PATTERN;idx=idx+1)begin
		@(negedge clk)
		valid = 1;
		in1 = in1_data[idx];
		in2 = in2_data[idx];
		in3 = in3_data[idx];
		in4 = in4_data[idx];
	end
@(negedge clk)
in1 = 0;
in2 = 0;
in3 = 0;
in4 = 0;
#(80*cycle) $finish;
end

endmodule
