`define INPUT_LENGTH 17
`define ANGLE_LENGTH 16
`define PIPELINE 2 
`define STAGE 13
module QR_decomposition_v2(clk, rst_n,valid, in_1, in_2, in_3, in_4, out_1, out_2, out_3, out_4);
	parameter HOLD = 4; //PIPELINE+2
	input clk;
	input rst_n;
	input valid; //for counter to start count up.
	input signed [`INPUT_LENGTH-1:0] in_1, in_2, in_3, in_4;
	output signed [`INPUT_LENGTH-1:0] out_1, out_2, out_3, out_4;

	// for clock counter and mode_in
	reg [2:0] counter;
	always @(posedge clk, negedge rst_n)begin
		if(!rst_n)begin
			//counter <= -(HOLD%8-1);
			counter <= 0;
		end
		else begin
			//counter <= (valid)? (counter==(8+3*HOLD-1))? 6'd0 : counter+1'b1 : counter;
			counter <= (valid)? (counter==3'd7)? 6'd0 : counter+1'b1 : counter;
		end
	end

	// for input data skew
	reg signed [`INPUT_LENGTH-1:0] in_skew_3[1*HOLD-1:0], in_skew_4[2*HOLD-1:0];
	reg signed [`INPUT_LENGTH-1:0] DU2[HOLD-1:0], DU3[HOLD-1:0];
	
	// data skew
	always @(posedge clk, negedge rst_n)begin
		if(!rst_n)begin
			in_skew_3[0] <= `INPUT_LENGTH'd0;
			in_skew_4[0] <= `INPUT_LENGTH'd0;
		end
		else begin
			in_skew_3[0] <= in_3;
			in_skew_4[0] <= in_4;
		end
	end
	generate
		genvar idx;
		for(idx=1; idx<1*HOLD; idx=idx+1)begin : in_skew_3_save
			always @(posedge clk, negedge rst_n)begin
				if(!rst_n)begin
					in_skew_3[idx] <= `INPUT_LENGTH'd0;
				end
				else begin
					in_skew_3[idx] <= in_skew_3[idx-1];
				end
			end
		end
		for(idx=1; idx<2*HOLD; idx=idx+1)begin : in_skew_4_save
			always @(posedge clk, negedge rst_n)begin
				if(!rst_n)begin
					in_skew_4[idx] <= `INPUT_LENGTH'd0;
				end
				else begin
					in_skew_4[idx] <= in_skew_4[idx-1];
				end
			end
		end
	endgenerate

	// for first row of PE
	wire signed [`INPUT_LENGTH-1:0] xin_PE[6:1], yin_PE[6:1];
	wire signed [`INPUT_LENGTH-1:0] xout_PE[6:1], yout_PE[6:1];
	wire signed [`ANGLE_LENGTH-1:0] angle_in_PE[6:1], angle_out_PE[6:1];
	reg signed [`ANGLE_LENGTH-1:0] angle_dff_PE[6:1];
	wire map_flag_in_PE[6:1], map_flag_out_PE[6:1];
	reg map_flag_dff_PE[6:1];
	wire mode_in_PE[6:1], mode_out_PE[6:1];
	wire rote_zero_PE[6:1];
	assign xin_PE[1] = in_1;
	assign yin_PE[1] = in_2;
	assign xin_PE[2] = xout_PE[1];
	assign yin_PE[2] = in_skew_3[1*HOLD-1];
	assign xin_PE[3] = xout_PE[2];
	assign yin_PE[3] = in_skew_4[2*HOLD-1];
	assign mode_in_PE[1] = (counter ==         1'b0)? 1'b0: 1'b1;
	assign mode_in_PE[2] = (counter == ((1*HOLD)%8))? 1'b0: 1'b1;
	assign mode_in_PE[3] = (counter == ((2*HOLD)%8))? 1'b0: 1'b1;
	//assign rote_zero_PE[1] = 1'b0;
	//assign rote_zero_PE[2] = 1'b0;
	//assign rote_zero_PE[3] = 1'b0;

	generate
		for(idx=1;idx<=3;idx=idx+1)begin : first_row_PE
			CORDIC_PE first_row_PE (.clk(clk), .rst_n(rst_n), .mode(mode_in_PE[idx]), .x(xin_PE[idx]), .y(yin_PE[idx]),
				 .angle(angle_in_PE[idx]), .map_flag(map_flag_in_PE[idx]), .rote_zero_in(1'b0),
				 .x_out(xout_PE[idx]), .y_out(yout_PE[idx]), 
				 .angle_out(angle_out_PE[idx]), .map_flag_out(map_flag_out_PE[idx]), .mode_out(mode_out_PE[idx]));
		end
		for (idx=1;idx<=3;idx=idx+1)begin : first_row_angle
			assign angle_in_PE[idx] = (mode_in_PE[idx])? angle_dff_PE[idx] : $signed(0);
		end
		for (idx=1;idx<=3;idx=idx+1)begin : first_row_flag
			assign map_flag_in_PE[idx] = (mode_in_PE[idx])? map_flag_dff_PE[idx] : $signed(0);
		end
		for (idx=1;idx<=3;idx=idx+1)begin : first_row_angle_dff
			always@(posedge clk, negedge rst_n)begin
				if(!rst_n) begin
					angle_dff_PE[idx] <= `ANGLE_LENGTH'd0;
				end
				else begin
					angle_dff_PE[idx] <= (mode_out_PE[idx]==1'b1)? angle_dff_PE[idx] : angle_out_PE[idx];
				end
			end
		end
		for (idx=1;idx<=3;idx=idx+1)begin : first_row_flag_dff
			always@(posedge clk, negedge rst_n)begin
				if(!rst_n) begin
				 map_flag_dff_PE[idx] <= `ANGLE_LENGTH'd0;
				end
				else begin
				 map_flag_dff_PE[idx] <= (mode_out_PE[idx]==1'b1)? map_flag_dff_PE[idx] : map_flag_out_PE[idx];
				end
			end
		end
	endgenerate

	//  for second row and RU
	assign xin_PE[4] = DU2[HOLD-1];
	assign yin_PE[4] = yout_PE[2];
	assign xin_PE[5] = xout_PE[4];
	assign yin_PE[5] = yout_PE[3];
	assign mode_in_PE[4] = (counter == (((2*HOLD+1)%8)))? 1'b0: 1'b1;
	assign mode_in_PE[5] = (counter == (((3*HOLD+1)%8)))? 1'b0: 1'b1;
	assign rote_zero_PE[4] = (counter == ((2*HOLD)%8))? 1'b1: 1'b0;
	assign rote_zero_PE[5] = (counter == ((3*HOLD)%8))? 1'b1: 1'b0;

	always @(posedge clk, negedge rst_n)begin
		if(!rst_n)begin
			DU2[0] <= `INPUT_LENGTH'd0;
		end
		else begin
			DU2[0] <= yout_PE[1];
		end
	end
	generate 
		for(idx=1; idx<HOLD; idx=idx+1)begin : DU2_save
			always @(posedge clk, negedge rst_n)begin
				if(!rst_n)begin
					DU2[idx] <= `INPUT_LENGTH'd0;
				end
				else begin
					DU2[idx] <= DU2[idx-1];
				end
			end
		end
		for(idx=4;idx<=5;idx=idx+1)begin : second_row_PE
			CORDIC_PE second_row_PE (.clk(clk), .rst_n(rst_n), .mode(mode_in_PE[idx]), .x(xin_PE[idx]), .y(yin_PE[idx]),
				 .angle(angle_in_PE[idx]), .map_flag(map_flag_in_PE[idx]), .rote_zero_in(rote_zero_PE[idx]),
				 .x_out(xout_PE[idx]), .y_out(yout_PE[idx]), 
				 .angle_out(angle_out_PE[idx]), .map_flag_out(map_flag_out_PE[idx]), .mode_out(mode_out_PE[idx]));
		end
		for (idx=4;idx<=5;idx=idx+1)begin : second_row_angle
			assign angle_in_PE[idx] = (mode_in_PE[idx])? angle_dff_PE[idx] : $signed(0);
		end
		for (idx=4;idx<=5;idx=idx+1)begin : second_row_flag
			assign map_flag_in_PE[idx] = (mode_in_PE[idx])? map_flag_dff_PE[idx] : $signed(0);
		end
		for (idx=4;idx<=5;idx=idx+1)begin : second_row_angle_dff
			always@(posedge clk, negedge rst_n)begin
				if(!rst_n) begin
					angle_dff_PE[idx] <= `ANGLE_LENGTH'd0;
				end
				else begin
					angle_dff_PE[idx] <= (mode_out_PE[idx]==1'b1)? angle_dff_PE[idx] : angle_out_PE[idx];
				end
			end
		end
		for (idx=4;idx<=5;idx=idx+1)begin : second_row_flag_dff
			always@(posedge clk, negedge rst_n)begin
				if(!rst_n) begin
				 map_flag_dff_PE[idx] <= `ANGLE_LENGTH'd0;
				end
				else begin
				 map_flag_dff_PE[idx] <= (mode_out_PE[idx]==1'b1)? map_flag_dff_PE[idx] : map_flag_out_PE[idx];
				end
			end
		end
	endgenerate

	//  for third row and RU
	assign xin_PE[6] = DU3[HOLD-1];
	assign yin_PE[6] = yout_PE[5];
	assign mode_in_PE[6] =  (counter == ((4*HOLD+2)%8))? 1'b0: 1'b1;
	assign rote_zero_PE[6] = (counter ==  ((4*HOLD+1)%8) || counter == ((4*HOLD)%8))? 1'b1: 1'b0;

	always @(posedge clk, negedge rst_n)begin
		if(!rst_n)begin
			DU3[0] <= `INPUT_LENGTH'd0;
		end
		else begin
			DU3[0] <= yout_PE[4];
		end
	end
	generate 
		for(idx=1; idx<HOLD; idx=idx+1)begin : DU3_save
			always @(posedge clk, negedge rst_n)begin
				if(!rst_n)begin
					DU3[idx] <= `INPUT_LENGTH'd0;
				end
				else begin
					DU3[idx] <= DU3[idx-1];
				end
			end
		end
		for(idx=6;idx<=6;idx=idx+1)begin : third_row_PE
			CORDIC_PE third_row_PE (.clk(clk), .rst_n(rst_n), .mode(mode_in_PE[idx]), .x(xin_PE[idx]), .y(yin_PE[idx]),
				 .angle(angle_in_PE[idx]), .map_flag(map_flag_in_PE[idx]), .rote_zero_in(rote_zero_PE[idx]),
				 .x_out(xout_PE[idx]), .y_out(yout_PE[idx]), 
				 .angle_out(angle_out_PE[idx]), .map_flag_out(map_flag_out_PE[idx]), .mode_out(mode_out_PE[idx]));
		end
		for (idx=6;idx<=6;idx=idx+1)begin : third_row_angle
			assign angle_in_PE[idx] = (mode_in_PE[idx])? angle_dff_PE[idx] : $signed(0);
		end
		for (idx=6;idx<=6;idx=idx+1)begin : third_row_flag
			assign map_flag_in_PE[idx] = (mode_in_PE[idx])? map_flag_dff_PE[idx] : $signed(0);
		end
		for (idx=6;idx<=6;idx=idx+1)begin : third_row_angle_dff
			always@(posedge clk, negedge rst_n)begin
				if(!rst_n) begin
					angle_dff_PE[idx] <= `ANGLE_LENGTH'd0;
				end
				else begin
					angle_dff_PE[idx] <= (mode_out_PE[idx]==1'b1)? angle_dff_PE[idx] : angle_out_PE[idx];
				end
			end
		end
		for (idx=6;idx<=6;idx=idx+1)begin : third_row_flag_dff
			always@(posedge clk, negedge rst_n)begin
				if(!rst_n) begin
				 map_flag_dff_PE[idx] <= `ANGLE_LENGTH'd0;
				end
				else begin
				 map_flag_dff_PE[idx] <= (mode_out_PE[idx]==1'b1)? map_flag_dff_PE[idx] : map_flag_out_PE[idx];
				end
			end
		end
	endgenerate


	// final data skew and output
	reg signed [`INPUT_LENGTH-1:0] out_skew_2[1*HOLD-1:0], out_skew_1[2*HOLD-1:0];
	always @(posedge clk, negedge rst_n)begin
		if(!rst_n)begin
			out_skew_1[0] <= `INPUT_LENGTH'd0;
			out_skew_2[0] <= `INPUT_LENGTH'd0;
		end
		else begin
			out_skew_1[0] <= xout_PE[3];
			out_skew_2[0] <= xout_PE[5];
		end
	end	
	generate 
		for(idx=1; idx<1*HOLD; idx=idx+1)begin : out_skew_2_save
			always @(posedge clk, negedge rst_n)begin
				if(!rst_n)begin
					out_skew_2[idx] <= `INPUT_LENGTH'd0;
				end
				else begin
					out_skew_2[idx] <= out_skew_2[idx-1];
				end
			end
		end
		for(idx=1; idx<2*HOLD; idx=idx+1)begin : out_skew_1_save
			always @(posedge clk, negedge rst_n)begin
				if(!rst_n)begin
					out_skew_1[idx] <= `INPUT_LENGTH'd0;
				end
				else begin
					out_skew_1[idx] <= out_skew_1[idx-1];
				end
			end
		end
	endgenerate
		
	// output assign
	assign out_1 = out_skew_1[2*HOLD-1];
	assign out_2 = out_skew_2[1*HOLD-1];
	assign out_3 = xout_PE[6]; 
	assign out_4 = yout_PE[6];


endmodule
