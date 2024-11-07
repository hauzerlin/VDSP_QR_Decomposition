`ifndef INPUT_LENGTH
`define INPUT_LENGTH 17 
`endif

`ifndef ANGLE_LENGTH
`define ANGLE_LENGTH 16 
`endif

`ifndef PIPELINE
`define PIPELINE  2 
`endif

`ifndef STAGE
`define STAGE 13
`endif
// mode 0: Vectoring mode, mode 1: Rotation mode
module CORDIC_PE(clk, rst_n, mode, x, y, angle, map_flag, rote_zero_in, x_out, y_out, angle_out, map_flag_out, mode_out);
	// remember to change with `define parameter PIPELINE and STAGE
	parameter DISTANCE = 7; // *****$ceil(`STAGE/`PIPELINE);***** 
	parameter DIVIDE = 6; 	// *****STAGE % DISTANCE;*****

	input clk, rst_n;
	input mode, map_flag, rote_zero_in;
	input signed [`INPUT_LENGTH-1:0] x, y;
	input signed [`ANGLE_LENGTH-1:0] angle;
	output reg signed [`INPUT_LENGTH-1:0] x_out, y_out;
	output reg signed [`ANGLE_LENGTH-1:0] angle_out;
	output reg map_flag_out;
	output reg mode_out;
	
	// wire and regs
	// wire for accumlator input and output
		wire mode_in [`STAGE-1:0];
		wire map_in [`STAGE-1:0];
		wire rote_zero [`STAGE-1:0];
		reg mu_data[`STAGE-1:0];
		wire mu_in [`STAGE-1:0], mu_out[`STAGE-1:0];
		wire signed [`INPUT_LENGTH-1:0] x_acc_in  [`STAGE-1:0], y_acc_in  [`STAGE-1:0];
		wire signed [`INPUT_LENGTH-1:0] x_acc_out [`STAGE-1:0], y_acc_out [`STAGE-1:0];
		wire signed [`ANGLE_LENGTH-1:0] ang_acc_in  [`STAGE-1:0];
		wire signed [`ANGLE_LENGTH-1:0] ang_acc_out [`STAGE-1:0];
	// reg for pipeline data to save
		reg mode_data [`PIPELINE:0];
		reg map_data [`PIPELINE:0];
		reg rote_zero_data[`PIPELINE:0];
		reg signed [`INPUT_LENGTH-1:0] x_data [`PIPELINE:0], y_data [`PIPELINE:0];
		reg signed [`ANGLE_LENGTH-1:0] angle_data [`PIPELINE:0];

	// Initial Stage
	always @(posedge clk, negedge rst_n)begin
		if(!rst_n) begin
			mode_data[0] <= 1'b0;
			map_data[0] <= 1'b0;
			rote_zero_data[0] <= 1'b0;
			x_data[0] <= `INPUT_LENGTH'd0;
			y_data[0] <= `INPUT_LENGTH'd0;
			angle_data[0] <= `ANGLE_LENGTH'd0;
		end
		else begin
			mode_data[0] <= mode;
			map_data[0] <= (mode==1'b0)? (x<0)? 1'b1: 1'b0: map_data[0];
			rote_zero_data[0] <= rote_zero_in;
			x_data[0] <= (!rote_zero_in)? (mode==1'b0)? ((x<0)? -x: x): ((map_data[0]==1'b1)? -x: x): x;
			y_data[0] <= (!rote_zero_in)? (mode==1'b0)? y : ((map_data[0]==1'b1)? -y: y): y;
			angle_data[0] <= (!rote_zero_in)? (mode==1'b0)? `ANGLE_LENGTH'd0 : angle : `ANGLE_LENGTH'd0;
		end
	end
	// CSD and Output stage
	always @(posedge clk, negedge rst_n) begin
		if(!rst_n) begin
			mode_data[`PIPELINE] <= 1'b0;
			map_data[`PIPELINE] <= 1'b0;
			rote_zero_data[`PIPELINE] <= 1'b0;
			x_data[`PIPELINE] <= `INPUT_LENGTH'd0;
			y_data[`PIPELINE] <= `INPUT_LENGTH'd0;
			angle_data[`PIPELINE] <= `ANGLE_LENGTH'd0;
		end
		else begin
			mode_data[`PIPELINE] <= mode_data[`PIPELINE-1];
			map_data[`PIPELINE] <= map_data[`PIPELINE-1];
			rote_zero_data[`PIPELINE] <= rote_zero_data[`PIPELINE-1];
			x_data[`PIPELINE] <= x_acc_out [`STAGE-1];
			y_data[`PIPELINE] <= y_acc_out [`STAGE-1];
			angle_data[`PIPELINE] <= ang_acc_out [`STAGE-1];
		end
	end
	wire signed [`INPUT_LENGTH-1:0] origin_x, origin_y, mapped_x, mapped_y;
	//assign origin_x = (!map_data[`PIPELINE])? -x_data[`PIPELINE]: x_data[`PIPELINE];
	//assign origin_y = (!map_data[`PIPELINE])? -y_data[`PIPELINE]: y_data[`PIPELINE];
	assign origin_x = x_data[`PIPELINE];
	assign origin_y = y_data[`PIPELINE];
	CSD scaling_x(.in(origin_x), .out(mapped_x));
	CSD scaling_y(.in(origin_y), .out(mapped_y));

	always @(posedge clk, negedge rst_n)begin
		if(!rst_n)begin
			x_out <= `INPUT_LENGTH'd0;
			y_out <= `INPUT_LENGTH'd0;
			angle_out <= `ANGLE_LENGTH'd0;
			map_flag_out <= 1'b0;
			mode_out <= 1'b1;
		end
		else begin
			x_out <= mapped_x; 
			y_out <= mapped_y;
			angle_out <= (mode_data[`PIPELINE]==1'b0)? (map_data[`PIPELINE])? -angle_data[`PIPELINE] : angle_data[`PIPELINE] : angle_data[`PIPELINE] ;
			map_flag_out <= map_data[`PIPELINE];
			mode_out <= mode_data[`PIPELINE];
		end
	end

	// Caculate Stages
	assign mode_in[0] = mode_data[0];
	assign map_in[0] = map_data[0];
	assign rote_zero[0] = rote_zero_data[0];
	assign x_acc_in[0] = x_data[0];
	assign y_acc_in[0] = y_data[0];
	assign ang_acc_in[0] = angle_data[0];
	// define parameters
	parameter [`ANGLE_LENGTH*`STAGE-1:0] elementary_angle = {
	16'b0011001001000011,
	16'b0001110110101100,
	16'b0000111110101101,
	16'b0000011111110101,
	16'b0000001111111110,
	16'b0000000111111111,
	16'b0000000011111111,
	16'b0000000001111111,
	16'b0000000000111111,
	16'b0000000000011111,
	16'b0000000000001111,
	16'b0000000000000111,
	16'b0000000000000011};
		
	// generate Stages
		genvar idx;
		generate
			for(idx=0; idx<`STAGE; idx=idx+1)begin
				CORDIC_PE_Acc #(.ELEMENTARY_ANGLE(elementary_angle[`ANGLE_LENGTH*(`STAGE-idx-1) +: `ANGLE_LENGTH]), .SHIFT_STAGE(idx)) 
				  Stage_(.x(x_acc_in[idx]), .y(y_acc_in[idx]), .ang(ang_acc_in[idx]), .mode(mode_in[idx]), .mu_in(mu_in[idx]), .map_in(map_in[idx]), .rote_zero(rote_zero[idx]),
				         .x_out(x_acc_out[idx]), .y_out(y_acc_out[idx]),.ang_out(ang_acc_out[idx]),.mu_out(mu_out[idx]));
			end
			for(idx=0;idx<`STAGE;idx=idx+1)begin : mu_data_save
				always @(posedge clk, negedge rst_n)begin
					if(!rst_n) begin
						mu_data[idx] <= 1'b0;
					end
					else begin
						mu_data[idx] <= mu_out[idx];
					end
				end
				assign mu_in[idx] = mu_data[idx];
			end
			for(idx=1; idx<`STAGE; idx=idx+1)begin
				if(((idx%DISTANCE) == 0)) begin
					assign mode_in[idx] = mode_data[(idx/DISTANCE)];
					assign map_in[idx] = map_data[(idx/DISTANCE)];
					assign rote_zero[idx] = rote_zero_data[(idx/DISTANCE)];
					assign x_acc_in[idx] = x_data[(idx/DISTANCE)];
					assign y_acc_in[idx] = y_data[(idx/DISTANCE)];
					assign ang_acc_in[idx] = angle_data[(idx/DISTANCE)];
				end
				else begin
					assign mode_in[idx] = mode_data[((idx-(idx%DISTANCE))/DISTANCE)];
					assign map_in[idx] = map_data[((idx-(idx%DISTANCE))/DISTANCE)];
					assign rote_zero[idx] = rote_zero_data[((idx-(idx%DISTANCE))/DISTANCE)];
					assign x_acc_in[idx] = x_acc_out[idx-1];
					assign y_acc_in[idx] = y_acc_out[idx-1];
					assign ang_acc_in[idx] = ang_acc_out[idx-1];
				end
			end
			for(idx=0; idx<`STAGE; idx=idx+1)begin
				if(!((idx==`STAGE-1)&&(DIVIDE==0)))begin
					if((idx)%DISTANCE== (DISTANCE-1)) begin
						always @(posedge clk, negedge rst_n)begin
							if(!rst_n)begin
								mode_data[(idx+1)/DISTANCE] <= 1'b0;
								map_data[(idx+1)/DISTANCE] <= 1'b0;
								rote_zero_data[(idx+1)/DISTANCE] <= 1'b0;
								x_data[(idx+1)/DISTANCE] <= `INPUT_LENGTH'd0;		
								y_data[(idx+1)/DISTANCE] <= `INPUT_LENGTH'd0;		
								angle_data[(idx+1)/DISTANCE] <= `ANGLE_LENGTH'd0;		
							end
							else begin
								mode_data[(idx+1)/DISTANCE] <= mode_data[(idx+1)/DISTANCE-1];
								map_data[(idx+1)/DISTANCE] <= map_data[(idx+1)/DISTANCE-1];
								rote_zero_data[(idx+1)/DISTANCE] <= rote_zero_data[(idx+1)/DISTANCE-1];
								x_data[(idx+1)/DISTANCE] <= x_acc_out[(idx)];
								y_data[(idx+1)/DISTANCE] <= y_acc_out[(idx)];
								angle_data[(idx+1)/DISTANCE] <= ang_acc_out[(idx)];
							end
						end
					end
				end
			end
		endgenerate
endmodule
