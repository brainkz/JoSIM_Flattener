



from itertools import product


period = list(range(60, 180, 10))
bias_scale = list(range(60, 180, 10))


    header = '''
.include /Users/brainkz/Documents/GitHub/SFQ_gate_Exploration/jj_library.cir
.include /Users/brainkz/Documents/GitHub/SFQ/sunmagnetics_library.cir
'''



    circuit = '''
.subckt	pg_init	a	b	clk	pout	gout
  * Need	3	clocks	here
  X_spl_a	a	a1	a2	spl2 BiasScale=0.9*dcbias
  X_spl_b	b	b1	b2	spl2 BiasScale=0.9*dcbias

  X_dff_a	a1	clk	a1_sync	dff BiasScale=dcbias BiasScale=0.9*dcbias
  X_dff_b	b1	clk	b1_sync	dff BiasScale=dcbias BiasScale=0.9*dcbias
  
  X_xor	a2	b2	clk	pout	xor BiasScale=dcbias
  X_and	a1_sync	b1_sync	gout	merge	mbias=0.2 BiasScale=dcbias
.ends

.subckt	B_pg	P0  G0  P1  G1  clk pout    gout
  * Need 4 clocks here
  X_spl_g1	G1	G1_copy_1	G1_copy_2	spl2 BiasScale=0.9*dcbias
  X_spl_p1	P1	P1_copy_1	P1_copy_2	spl2 BiasScale=0.9*dcbias
  
  X_pg	P1_copy_1	G1_copy_1	pg	merge	mbias=1 BiasScale=dcbias
  X_gg	G0          G1_copy_2	gg	merge	mbias=1 BiasScale=dcbias

  X_dff_p0	P0	        clk	P0_sync	dff BiasScale=dcbias
  X_dff_p1	P1_copy_2	clk	P1_sync	dff BiasScale=dcbias
  X_dff_pg    pg          clk pg_sync dff BiasScale=dcbias
  X_dff_gg    gg          clk gg_sync dff BiasScale=dcbias

  X_and_g	pg_sync	gg_sync	gout	merge	mbias=0.2 BiasScale=dcbias
  X_and_p	P0_sync	P1_sync	pout	merge	mbias=0.2 BiasScale=dcbias
.ends

.subckt	W_pg	pin  gin  clk pout    gout
  * Need 2 clocks here
  Xp	    pin     clk	    pout	dff BiasScale=dcbias
  Xg	    gin     clk	    gout	dff BiasScale=dcbias
.ends

.subckt	B_g         G0  P1  G1  clk          gout
  * Need 2 clocks here
  X_spl_g1	G1	G1_copy_1	G1_copy_2	spl2 BiasScale=0.9*dcbias
  
  X_pg	P1   G1_copy_1	pg	merge	mbias=1 BiasScale=dcbias
  X_gg	G0   G1_copy_2	gg	merge	mbias=1 BiasScale=dcbias

  X_dff_pg    pg          clk pg_sync dff BiasScale=dcbias
  X_dff_gg    gg          clk gg_sync dff BiasScale=dcbias

  X_and_g	pg_sync	gg_sync	gout	merge	mbias=0.2 BiasScale=dcbias
.ends

.subckt ptl a q
    X_tx  a       a_ptl LSmitll_PTLTX BiasScale=dcbias
    X_rx  a_ptl   q     LSmitll_PTLRX BiasScale=dcbias
.ends

.param tclock=180e-12
.param OS=tclock/40
.param step=0.08
*  Input signals
  Xa0 a0 data_sig_a clk_offset=OS do_factor=step*2+step*0 Tck=tclock
  Xb0 b0 data_sig_b clk_offset=OS do_factor=step*3+step*0 Tck=tclock
  Xa1 a1 data_sig_c clk_offset=OS do_factor=step*4+step*0 Tck=tclock
  Xb1 b1 data_sig_d clk_offset=OS do_factor=step*5+step*0 Tck=tclock
  Xa2 a2 data_sig_e clk_offset=OS do_factor=step*6+step*0 Tck=tclock
  Xb2 b2 data_sig_f clk_offset=OS do_factor=step*7+step*0 Tck=tclock
  Xa3 a3 data_sig_g clk_offset=OS do_factor=step*8+step*0 Tck=tclock
  Xb3 b3 data_sig_h clk_offset=OS do_factor=step*9+step*0 Tck=tclock
* 
* Cycle 0
  xt00 T00 clk_sig clk_offset=OS*4 factor=3 Tck=tclock
  Xi0 a0 b0 T00 iP0_0 iG0_0 pg_init
  
  Xt01 T01 clk_sig clk_offset=OS*4 factor=3 Tck=tclock
  Xi1 a1 b1 T01 iP1_0 iG1_0 pg_init
  
  Xt02 T02 clk_sig clk_offset=OS*4 factor=3 Tck=tclock
  Xi2 a2 b2 T02 iP2_0 iG2_0 pg_init
  
  Xt03 T03 clk_sig clk_offset=OS*4 factor=3 Tck=tclock
  Xi3 a3 b3 T03 iP3_0 iG3_0 pg_init
* 

*   PTL cycle 0
    X_ptl_iP0_0  iP0_0       iP0_0_rx   ptl
    X_ptl_iG0_0  iG0_0       iG0_0_rx   ptl
    X_ptl_iP1_0  iP1_0       iP1_0_rx   ptl
    X_ptl_iG1_0  iG1_0       iG1_0_rx   ptl
    X_ptl_iP2_0  iP2_0       iP2_0_rx   ptl
    X_ptl_iG2_0  iG2_0       iG2_0_rx   ptl
    X_ptl_iP3_0  iP3_0       iP3_0_rx   ptl
    X_ptl_iG3_0  iG3_0       iG3_0_rx   ptl
* 

* Splitters at level 0
  * Splitting 0
*   Xspl_iP0_0 iP0_0 iP0_0_to0 iP0_0_to1           spl2 BiasScale=0.9*dcbias
  Xspl_iG0_0 iG0_0_rx iG0_0_to0 iG0_0_to1           spl2 BiasScale=0.9*dcbias
  * Splitting 1
  Xspl_iP1_0 iP1_0_rx iP1_0_to1           iP1_0_out spl2 BiasScale=0.9*dcbias
  * Splitting 2
  Xspl_iP2_0 iP2_0_rx iP2_0_to2 iP2_0_to3 iP2_0_out spl3 BiasScale=0.8
  Xspl_iG2_0 iG2_0_rx iG2_0_to2 iG2_0_to3           spl2 BiasScale=0.9*dcbias
  * Splitting 3
  Xspl_iP3_0 iP3_0_rx iP3_0_to1           iP3_0_out spl2 BiasScale=0.9*dcbias
* 

*   
  *   R_ip0_to0     iP0_0_to0  0   1
  *   R_ip0_to1     iP0_0_to1  0   1
*   R_ig0_to0     iG0_0_to0  0   1
*   R_ig0_to1     iG0_0_to1  0   1
*   R_ip1_to1     iP1_0_to1  0   1
*   R_ip1_out     iP1_0_out  0   1
*   R_ig1_rx      iG1_0_rx   0   1
*   R_ip2_to2     iP2_0_to2  0   1
*   R_ip2_to3     iP2_0_to3  0   1
*   R_ip2_out     iP2_0_out  0   1
*   R_ig2_to2     iG2_0_to2  0   1
*   R_ig2_to3     iG2_0_to3  0   1
*   R_ip3_to1     iP3_0_to1  0   1
*   R_ip3_out     iP3_0_out  0   1
*   R_ig3         iG3_0_rx  0   1

*  Cycle 1
  Xt04 T04 clk_sig clk_offset=OS*3 factor=2 Tck=tclock
*   X_pg0_01 iP0_0_to0 iG0_0_to0                 T04 P0_1 G0_1 W_pg
  X_pg0_01 iP0_0_rx     iG0_0_to0                    T04 P0_1 G0_1 W_pg

  Xt05 T05 clk_sig clk_offset=OS*3 factor=2 Tck=tclock
*   X_pg1_01 iP0_0_to1 iG0_0_to1 iP1_0_to1 iG1_0 T05 P1_1 G1_1 B_pg
  X_pg1_01              iG0_0_to1 iP1_0_to1 iG1_0_rx T05      G1_1 B_g

  Xt06 T06 clk_sig clk_offset=OS*3 factor=2 Tck=tclock
  X_pg2_01 iP2_0_to2    iG2_0_to2                    T06 P2_1 G2_1 W_pg

  Xt07 T07 clk_sig clk_offset=OS*3 factor=4 Tck=tclock
  X_pg3_01 iP2_0_to3    iG2_0_to3 iP3_0_to1 iG3_0_rx T07 P3_1 G3_1 B_pg

  Xd01 D01 clk_sig clk_offset=OS*3 factor=1 Tck=tclock
  X_dff_ip1_01  iP1_0_out   D01 iP1_1_out   dff BiasScale=dcbias

  Xd02 D02 clk_sig clk_offset=OS*3 factor=1 Tck=tclock
*   Xbuf_ip2_01   iP2_0_out       iP2_0_tmp     LSMITLL_BUFF
  X_dff_ip2_01  iP2_0_out   D02 iP2_1_out   dff BiasScale=dcbias

  Xd03 D03 clk_sig clk_offset=OS*3 factor=1 Tck=tclock
  X_dff_ip3_01  iP3_0_out   D03 iP3_1_out   dff BiasScale=dcbias
* 

*   PTL cycle 1
    X_ptl_P0_1  P0_1       P0_1_rx   ptl
    X_ptl_G0_1  G0_1       G0_1_rx   ptl

    X_ptl_G1_1  G1_1       G1_1_rx   ptl
    X_ptl_P2_1  P2_1       P2_1_rx   ptl
    X_ptl_G2_1  G2_1       G2_1_rx   ptl
    X_ptl_P3_1  P3_1       P3_1_rx   ptl
    X_ptl_G3_1  G3_1       G3_1_rx   ptl

    X_ptl_iP1_1  iP1_1_out       iP1_1_out_rx   ptl
    X_ptl_iP2_1  iP2_1_out       iP2_1_out_rx   ptl
    X_ptl_iP3_1  iP3_1_out       iP3_1_out_rx   ptl
* 
* Splitters at level 1
  * Splitting 0
  * Splitting 1
*   Xspl_P1_1   P1_1    P1_1_to1    P1_1_to2    P1_1_to3    spl3 BiasScale=0.8
  Xspl_G1_1   G1_1_rx    G1_1_to1    G1_1_to2    G1_1_to3    spl3 BiasScale=0.8
  * Splitting 2
  * Splitting 3
* 

*  Cycle 2
  Xt08 T08 clk_sig clk_offset=OS*2 factor=2 Tck=tclock
  X_pg0_12  P0_1_rx      G0_1_rx                 T08 P0_2 G0_2 W_pg

  Xt09 T09 clk_sig clk_offset=OS*2 factor=1 Tck=tclock
*   X_pg1_12  P1_1_to1  G1_1_to1             T09 P1_2 G1_2 W_pg
  X_pg1_12            G1_1_to1             T09      G1_2 dff BiasScale=dcbias

  Xt10 T10 clk_sig clk_offset=OS*2 factor=2 Tck=tclock
*   X_pg2_12  P1_1_to2  G1_1_to2  P2_1  G2_1 T10 P2_2 G2_2 B_pg
  X_pg2_12            G1_1_to2  P2_1_rx  G2_1_rx T10      G2_2 B_g

  Xt11 T11 clk_sig clk_offset=OS*2 factor=2 Tck=tclock
*   X_pg3_12  P1_1_to3  G1_1_to3  P3_1  G3_1 T11 P3_2 G3_2 B_pg
  X_pg3_12            G1_1_to3  P3_1_rx  G3_1_rx T11      G3_2 B_g

  Xd11 D11 clk_sig clk_offset=OS*2 factor=1 Tck=tclock
  X_dff_ip1_12  iP1_1_out_rx   D11 iP1_2_out   dff BiasScale=dcbias

  Xd12 D12 clk_sig clk_offset=OS*2 factor=1 Tck=tclock
*   Xbuf_ip2_12   iP2_1_out       iP2_1_tmp     LSMITLL_BUFF
  X_dff_ip2_12  iP2_1_out_rx   D12 iP2_2_out   dff BiasScale=dcbias

  Xd13 D13 clk_sig clk_offset=OS*2 factor=1 Tck=tclock
  X_dff_ip3_12  iP3_1_out_rx   D13 iP3_2_out   dff BiasScale=dcbias 
* 
*   PTL cycle 1
    X_ptl_P0_2  P0_2       P0_2_rx   ptl
    X_ptl_G0_2  G0_2       G0_2_rx   ptl

    X_ptl_G1_2  G1_2       G1_2_rx   ptl
    * X_ptl_P2_1  P2_1       P2_1_rx   ptl
    X_ptl_G2_2  G2_2       G2_2_rx   ptl
    * X_ptl_P3_1  P3_1       P3_1_rx   ptl
    X_ptl_G3_2  G3_2       G3_2_rx   ptl

    X_ptl_iP1_2  iP1_2_out       iP1_2_out_rx   ptl
    X_ptl_iP2_2  iP2_2_out       iP2_2_out_rx   ptl
    X_ptl_iP3_2  iP3_2_out       iP3_2_out_rx   ptl
* 

* Final processing
  * Bit 0
  Xt12 T12 clk_sig clk_offset=OS*1 factor=1 Tck=tclock
  X_s0  P0_2_rx                T12 s0  dff BiasScale=dcbias
  R_s0  s0      0   1
*   R_c1  G0_2    0   1
  Xt13 T13 clk_sig clk_offset=OS*1 factor=1 Tck=tclock
  Xbuf_G0_2_tmp   G0_2_rx         G0_2_tmp     LSMITLL_BUFF
  Xbuf_G0_2   G0_2_tmp     G0_2_buf     jtl
  Xbuf_iP1_2_tmp iP1_2_out_rx    iP1_2_tmp     LSMITLL_BUFF
  Xbuf_iP1_2 iP1_2_tmp    iP1_2_out_buf jtl
  X_s1	G0_2_buf	iP1_2_out_buf	T13	s1	xor t_scaling=1 main_scaling=0.8 bottom_scaling=0.8
  R_s1  s1      0   1
*   R_c2  G1_2    0   1
  Xt14 T14 clk_sig clk_offset=OS*1 factor=1 Tck=tclock
  Xbuf_G1_2_tmp   G1_2_rx         G1_2_tmp     LSMITLL_BUFF
  Xbuf_G1_2   G1_2_tmp     G1_2_buf     jtl
  Xbuf_iP2_2_tmp iP2_2_out_rx    iP2_2_tmp     LSMITLL_BUFF
  Xbuf_iP2_2 iP2_2_tmp    iP2_2_out_buf jtl
  X_s2	G1_2_buf	iP2_2_out_buf	T14	s2	xor t_scaling=1 main_scaling=0.8 bottom_scaling=0.8
  R_s2  s2    0   1
*   R_c3  G2_2    0   1
  Xt15 T15 clk_sig clk_offset=OS*1 factor=1 Tck=tclock
*   Xbuf_G2_2_tmp   G2_2_rx         G2_2_tmp     LSMITLL_BUFF
*   Xbuf_G2_2   G2_2_tmp     G2_2_buf     jtl
*   Xbuf_iP3_2_tmp iP3_2_out    iP3_2_tmp     LSMITLL_BUFF
*   Xbuf_iP3_2 iP3_2_tmp    iP3_2_out_buf jtl
  X_s3	G2_2_rx	iP3_2_out_rx	T15	s3	xor t_scaling=1 main_scaling=0.8 bottom_scaling=0.8
  R_s3  s3    0   1
*   R_c4  G3_2    0   1
  Xt16 T16 clk_sig clk_offset=OS*1 factor=1 Tck=tclock
  X_s4  G3_2_rx                T16 s4  dff BiasScale=dcbias
  R_s4  s4      0   1
* 
    * R_p0_2     P0_2  0   1
    * R_g0_2     G0_2  0   1
    * R_p1_2     P1_2  0   1
    * R_g1_2     G1_2  0   1
    * R_p2_2     P2_2  0   1
    * R_g2_2     G2_2  0   1
    * R_p3_2     P3_2  0   1
    * R_g3_2     G3_2  0   1
    * R_ip1_2     iP1_2_out_rx  0   1
    * R_ip2_2     iP2_2_out_rx  0   1
    * R_ip3_2     iP3_2_out_rx  0   1
    * .tran 1e-12 52000e-12
    * .tran 1e-12 39000e-12
    .tran 1e-12 52000e-12/200*180
    .end 


* *   
*   R_ip0      P0_1  0   1
*   R_ig0      G0_1  0   1
*   R_ip1     iP1_0  0   1
*   R_ig1     iG1_0  0   1
*   R_ip2     iP2_0  0   1
*   R_ig2     iG2_0  0   1
*   R_ip3     iP3_0  0   1
*   R_ig3     iG3_0  0   1
* * 
'''