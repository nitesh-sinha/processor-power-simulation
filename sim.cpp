#include "common.h"
#include "sim.h"
#include "trace.h" 
#include <stdlib.h>
#include "cache.h"  /**** NEW-LAB2*/ 
#include "memory.h" // NEW-LAB2 
#include <ctype.h> /* Library for useful character operations */
#include "bpred.h"
#include "vmem.h"
#include "knob.h"
#include "all_knobs.h"

/*******************************************************************/
/* Simulator frame */ 
/*******************************************************************/

bool run_a_cycle(memory_c *m ); // please modify run_a_cycle function argument  /** NEW-LAB2 */ 
void init_structures(memory_c *m); // please modify init_structures function argument  /** NEW-LAB2 */ 

/* uop_pool related variables */ 

uint32_t free_op_num;
uint32_t active_op_num; 
Op *op_pool; 
Op *op_pool_free_head = NULL; 
Op* dummy_op;
int ghr_len;
uint64_t cd_stall_addr;


bpred *branchpred;
tlb *dtlb;

// knob variables
KnobsContainer *g_knobsContainer; /* < knob container > */
all_knobs_c    *g_knobs; /* < all knob variables > */

gzFile g_stream;

/* simulator related local functions */ 

bool icache_access(ADDRINT addr);
bool dcache_access(ADDRINT addr);
void  init_latches(void);



void init_knobs(int argc, char** argv)
{
  // Create the knob managing class
  g_knobsContainer = new KnobsContainer();

  // Get a reference to the actual knobs for this component instance
  g_knobs = g_knobsContainer->getAllKnobs();

  // apply the supplied command line switches
  char* pInvalidArgument = NULL;
  g_knobsContainer->applyComandLineArguments(argc, argv, &pInvalidArgument);

  g_knobs->display();
}

void read_trace_file(void)
{
  g_stream = gzopen((KNOB(KNOB_TRACE_FILE)->getValue()).c_str(), "r");
}

// simulator main function is called from outside of this file 

void simulator_main(int argc, char** argv) 
{
  init_knobs(argc, argv);

  // trace driven simulation 
  read_trace_file();
   /** NEW-LAB2 */ /* just note: passing main memory pointers is hack to mix up c++ objects and c-style code */  /* Not recommended at all */ 
  memory_c *main_memory = new memory_c();  // /** NEW-LAB2 */ 
  init_structures(main_memory);  // please modify run_a_cycle function argument  /** NEW-LAB2 */ 
  run_a_cycle(main_memory); // please modify run_a_cycle function argument  /** NEW-LAB2 */ 

}
int op_latency[NUM_OP_TYPE]; 

void init_op_latency(void)
{
  op_latency[OP_INV]   = 1; 
  op_latency[OP_NOP]   = 1; 
  op_latency[OP_CF]    = 1; 
  op_latency[OP_CMOV]  = 1; 
  op_latency[OP_LDA]   = 1;
  op_latency[OP_LD]    = 1; 
  op_latency[OP_ST]    = 1; 
  op_latency[OP_IADD]  = 1; 
  op_latency[OP_IMUL]  = 2; 
  op_latency[OP_IDIV]  = 4; 
  op_latency[OP_ICMP]  = 2; 
  op_latency[OP_LOGIC] = 1; 
  op_latency[OP_SHIFT] = 2; 
  op_latency[OP_BYTE]  = 1; 
  op_latency[OP_MM]    = 2; 
  op_latency[OP_FMEM]  = 2; 
  op_latency[OP_FCF]   = 1; 
  op_latency[OP_FCVT]  = 4; 
  op_latency[OP_FADD]  = 2; 
  op_latency[OP_FMUL]  = 4; 
  op_latency[OP_FDIV]  = 16; 
  op_latency[OP_FCMP]  = 2; 
  op_latency[OP_FBIT]  = 2; 
  op_latency[OP_FCMP]  = 2; 
}

void init_op(Op *op)
{
  op->num_src               = 0; 
  op->src[0]                = -1; 
  op->src[1]                = -1;
  op->dst                   = -1; 
  op->opcode                = 0; 
  op->is_fp                 = false;
  op->cf_type               = NOT_CF;
  op->mem_type              = NOT_MEM;
  op->write_flag             = 0;
  op->inst_size             = 0;
  op->ld_vaddr              = 0;
  op->st_vaddr              = 0;
  op->instruction_addr      = 0;
  op->branch_target         = 0;
  op->actually_taken        = 0;
  op->mem_read_size         = 0;
  op->mem_write_size        = 0;
  op->valid                 = FALSE; 
  op->last_op		    	= FALSE; // initialize last_op to FALSE
  op->dummy					=false;
  /* you might add more features here */ 
}


void init_op_pool(void)
{
  /* initialize op pool */ 
  op_pool = new Op [1024];
  free_op_num = 1024; 
  active_op_num = 0; 
  uint32_t op_pool_entries = 0; 
  int ii;
  for (ii = 0; ii < 1023; ii++) {

    op_pool[ii].op_pool_next = &op_pool[ii+1]; 
    op_pool[ii].op_pool_id   = op_pool_entries++; 
    init_op(&op_pool[ii]); 
  }
  op_pool[ii].op_pool_next = op_pool_free_head; 
  op_pool[ii].op_pool_id   = op_pool_entries++;
  init_op(&op_pool[ii]); 
  op_pool_free_head = &op_pool[0]; 
}





Op *get_free_op(void)
{
  /* return a free op from op pool */ 

  if (op_pool_free_head == NULL || (free_op_num == 1)) {
    std::cout <<"ERROR! OP_POOL SIZE is too small!! " << endl; 
    std::cout <<"please check free_op function " << endl; 
    assert(1); 
    exit(1);
  }

  free_op_num--;
  assert(free_op_num); 

  Op *new_op = op_pool_free_head; 
  op_pool_free_head = new_op->op_pool_next; 
  assert(!new_op->valid); 
  init_op(new_op);
  active_op_num++; 
  return new_op; 
}

void free_op(Op *op)
{
  free_op_num++;
  active_op_num--; 
  op->valid = FALSE; 
  op->op_pool_next = op_pool_free_head;
  op_pool_free_head = op; 
}



/*******************************************************************/
/*  Data structure */
/*******************************************************************/

typedef struct pipeline_latch_struct {
  Op *op; /* you must update this data structure. */
  bool op_valid; // it tells a stage in the pipeline whether there is a valid instruction to execute in that cycle
   /* you might add more data structures. But you should complete the above data elements */ 
}pipeline_latch; 


typedef struct Reg_element_struct{
  bool valid; // 8 bit register
  // data is not needed 
  /* you might add more data structures. But you should complete the above data elements */ 
}REG_element; 

REG_element register_file[NUM_REG]; //NUM_REG is 32

// To initialize the 32 registers in the processor

void init_registers()
{
  for(int i=0;i<NUM_REG;i++)
    register_file[i].valid=TRUE;
}


/*******************************************************************/
/* These are the functions you'll have to write.  */ 
/*******************************************************************/

void FE_stage();
void ID_stage();
void EX_stage(); 
void MEM_stage(memory_c *main_memory); // please modify MEM_stage function argument  /** NEW-LAB2 */ 
void WB_stage(); 

/*******************************************************************/
/*  These are the variables you'll have to write.  */ 
/*******************************************************************/

bool sim_end_condition = FALSE;     /* please complete the condition. */ 
UINT64 retired_instruction = 0;    /* number of retired instruction. (only correct instructions) */ 
UINT64 cycle_count = 0;            /* total number of cycles */ 
UINT64 data_hazard_count = 0;  
UINT64 control_hazard_count = 0; 
UINT64 icache_miss_count = 0;      /* total number of icache misses. for Lab #2 and Lab #3 */ 
UINT64 dcache_hit_count = 0;      /* total number of dcache  misses. for Lab #2 and Lab #3 */ 
UINT64 dcache_miss_count = 0;      /* total number of dcache  misses. for Lab #2 and Lab #3 */ 
UINT64 l2_cache_miss_count = 0;    /* total number of L2 cache  misses. for Lab #2 and Lab #3 */  
UINT64 dram_row_buffer_hit_count = 0; /* total number of dram row buffer hit. for Lab #2 and Lab #3 */   // NEW-LAB2
UINT64 dram_row_buffer_miss_count = 0; /* total number of dram row buffer hit. for Lab #2 and Lab #3 */   // NEW-LAB2
UINT64 store_load_forwarding_count = 0;  /* total number of store load forwarding for Lab #2 and Lab #3 */  // NEW-LAB2

uint64_t bpred_mispred_count = 0;  /* total number of branch mispredictions */  // NEW-LAB3
uint64_t bpred_okpred_count = 0;   /* total number of correct branch predictions */  // NEW-LAB3
uint64_t dtlb_hit_count = 0;       /* total number of DTLB hits */ // NEW-LAB3
uint64_t dtlb_miss_count = 0;      /* total number of DTLB miss */ // NEW-LAB3


/* ********* Counters added for power simulation lab *************** */
uint64_t fp_instruction =0; /* Total number of Floating point instructions */
uint64_t int_instruction =0; /* Total number of integer instructions */
uint64_t branch_instruction =0; /* Total number of branch instructions */
uint64_t load_instruction =0; /* Total number of load instructions */
uint64_t store_instruction =0; /* Total number of store instructions */
uint64_t mul_instruction = 0; /* Total number of Multiply instructions */
uint64_t float_regfile_writes =0; /* Total number of FP register writes */
uint64_t int_regfile_writes =0; /* Total number of INT register writes */
uint64_t fp_sched_access = 0; /* Total number of scheduler access for FP instruction */
uint64_t int_sched_access = 0; /* Total number of scheduler access for Int instruction */
uint64_t fp_reg_read = 0; /* Total number of FP register reads */
uint64_t int_reg_read = 0; /* Total number of Int register reads */
uint64_t dcache_hit_load = 0; /* Total number of dcache hit for LOAD instruction */
uint64_t dcache_hit_store = 0; /* Total number of dcache hit for STORE instruction */
uint64_t inserted_in_dcache = 0; /* Total number of dcache inserts */



/* ********* Counters added for power simulation lab *************** */

pipeline_latch *MEM_latch;  
pipeline_latch *EX_latch;
pipeline_latch *ID_latch;
pipeline_latch *FE_latch;
UINT64 ld_st_buffer[LD_ST_BUFFER_SIZE]; 
UINT64 next_pc; 

Cache *data_cache = new Cache;  // NEW-LAB2 
int latency=0; // to keep track of instrn latency in EX stage

bool cd_stall=FALSE; // stall due to control dependency
bool dd_stall=FALSE; // stall due to data dependency
bool ex_stall=FALSE; // stall due to execution cycles
bool fe_stall=FALSE; // stall after last instrn fetched
bool mem_latency_stall=FALSE; // stall due to mem_latency > 1. for lab 2 
bool mshr_stall=FALSE; // for lab2
bool mshr1_stall=FALSE; // for lab3(MSHR stall for dummy op)
bool mshr2_stall=FALSE; // for lab3(MSHR stall for valid mem op)
bool mem_latency_flag=FALSE;
int dcache_latency;
bool finalEnd=false;

// flags to check if the current instruction in ID stage has control or data dependency
bool cd=FALSE;
bool dd=FALSE;


int fe_ctr=0; 				// to detect infinte loops in program

// To move instructions services by DRAM from op_queue to mem_op_queue for next cycle retirement
list <Op *> op_queue; 
list <Op *> mem_op_queue;

bool last_op_in_queue=false; // it will be true when last instruction reaches mem_op_queue
bool use_bpred, use_vmem; 	 // Whether branch prediction and virtual memory is used in pipeline or not
bool btaken; 				 // Branch predicted direction
bool tlb_stall=false; 		 // for blocking TLB implementation
bool dummy_op_flag=false; 	 // to catch the dummy op after it comes back from DRAM

bool tlb_inst=false;		 // to make sure tlb_install and tlb_access are not done in same cycle
uint64_t vpn,pfn,pte,paddr;
int vpn_page_size;

/*******************************************************************/
/*  Print messages  */
/*******************************************************************/

void print_stats() {
  std::ofstream out((KNOB(KNOB_OUTPUT_FILE)->getValue()).c_str());
  /* Do not modify this function. This messages will be used for grading */ 
  out << "Total instruction: " << retired_instruction << endl; 
  out << "Total cycles: " << cycle_count << endl; 
  float ipc = (cycle_count ? ((float)retired_instruction/(float)cycle_count): 0 );
  out << "IPC: " << ipc << endl; 
  out << "Total I-cache miss: " << icache_miss_count << endl; 
  out << "Total D-cache hit: " << dcache_hit_count << endl; 
  out << "Total D-cache miss: " << dcache_miss_count << endl; 
  out << "Total L2-cache miss: " << l2_cache_miss_count << endl; 
  out << "Total data hazard: " << data_hazard_count << endl;
  out << "Total control hazard : " << control_hazard_count << endl; 
  out << "Total DRAM ROW BUFFER Hit: " << dram_row_buffer_hit_count << endl;  // NEW-LAB2
  out << "Total DRAM ROW BUFFER Miss: "<< dram_row_buffer_miss_count << endl;  // NEW-LAB2 
  out <<" Total Store-load forwarding: " << store_load_forwarding_count << endl;  // NEW-LAB2 

  // new for LAB3
  out << "Total Branch Predictor Mispredictions: " << bpred_mispred_count << endl;
  out << "Total Branch Predictor OK predictions: " << bpred_okpred_count << endl;
  out << "Total DTLB Hit: " << dtlb_hit_count << endl;
  out << "Total DTLB Miss: " << dtlb_miss_count << endl;

  // new for LAB4(17 new counters)

  out << "Total number of Floating point instructions: " << fp_instruction << endl;
  out << "Total number of integer instructions: " << int_instruction << endl;
  out << "Total number of branch instructions: " << branch_instruction << endl;
  out << "Total number of load instructions: " << load_instruction << endl;
  out << "Total number of store instructions: " << store_instruction << endl;
  out << "Total number of Multiply instructions: " << mul_instruction << endl;
  out << "Total number of FP register writes: " << float_regfile_writes << endl;
  out << "Total number of INT register writes: " << int_regfile_writes << endl;
  out << "Total number of scheduler access for FP instruction: " << fp_sched_access << endl;
  out << "Total number of scheduler access for Int instruction: " << int_sched_access << endl;
  out << "Total number of ROB reads/writes: " << int_sched_access + fp_sched_access << endl;
  out << "Total number of register reads for Int ops: " << int_reg_read << endl;
  out << "Total number of register reads for FP ops: " << fp_reg_read << endl;
  out << "Total DTLB accesses(read+write): " << (dtlb_hit_count + 2*dtlb_miss_count) << endl;
  out << "Total dcache reads: " << (dcache_hit_count - dcache_hit_store) << endl;
  out << "Total dcache writes: " << (inserted_in_dcache + dcache_hit_store) << endl;
  out << "Total memory accesses: " << (dram_row_buffer_hit_count + dram_row_buffer_miss_count) << endl;

  // check for dcache writes


  out.close();
}


/*******************************************************************/
/*  Support Functions  */ 
/*******************************************************************/

bool get_op(Op *op)
{
  static UINT64 unique_count = 0; 
  Trace_op trace_op; 
  bool success = FALSE; 
  // read trace 
  // fill out op info 
  // return FALSE if the end of trace 
  success = (gzread(g_stream, &trace_op, sizeof(Trace_op)) == sizeof(Trace_op));
  if (KNOB(KNOB_PRINT_INST)->getValue()) dprint_trace(&trace_op); 

  /* copy trace structure to op */ 
  if (success) { 
    copy_trace_op(&trace_op, op); 

    op->inst_id  = unique_count++;
    op->valid    = TRUE; 
  }
  return success; 
}
/* return op execution cycle latency */ 

int get_op_latency(Op *op) 
{
  assert (op->opcode < NUM_OP_TYPE); 
  return op_latency[op->opcode];
}

/* Print out all the register values */ 
void dump_reg() {
  for (int ii = 0; ii < NUM_REG; ii++) {
    std::cout << cycle_count << ":register[" << ii  << "]: V:" << register_file[ii].valid << endl; 
  }
}

void print_pipeline() {
  std::cout << "--------------------------------------------" << endl; 
  std::cout <<"cycle count : " << dec << cycle_count << " retired_instruction : " << retired_instruction << endl; 
  std::cout << (int)cycle_count << " FE: " ;
  if (FE_latch->op_valid) {
    Op *op = FE_latch->op; 
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }
  std::cout << " ID: " ;
  if (ID_latch->op_valid) {
    Op *op = ID_latch->op; 
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }
  std::cout << " EX: " ;
  if (EX_latch->op_valid) {
    Op *op = EX_latch->op; 
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }


  std::cout << " MEM: " ;
  if (MEM_latch->op_valid) {
    Op *op = MEM_latch->op; 
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }
  cout << endl; 
  //  dump_reg();   
  std::cout << "--------------------------------------------" << endl; 
}

void print_heartbeat()
{
  static uint64_t last_cycle ;
  static uint64_t last_inst_count; 
  float temp_ipc = float(retired_instruction - last_inst_count) /(float)(cycle_count-last_cycle) ;
  float ipc = float(retired_instruction) /(float)(cycle_count) ;
  /* Do not modify this function. This messages will be used for grading */ 
  cout <<"**Heartbeat** cycle_count: " << cycle_count << " inst:" << retired_instruction << " IPC: " << temp_ipc << " Overall IPC: " << ipc << endl; 
  last_cycle = cycle_count;
  last_inst_count = retired_instruction; 
}
/*******************************************************************/
/*                                                                 */
/*******************************************************************/

bool run_a_cycle(memory_c *main_memory){  // please modify run_a_cycle function argument  /** NEW-LAB2 */ 

  for (;;) { 
    if (((KNOB(KNOB_MAX_SIM_COUNT)->getValue() && (cycle_count >= KNOB(KNOB_MAX_SIM_COUNT)->getValue())) || 
      (KNOB(KNOB_MAX_INST_COUNT)->getValue() && (retired_instruction >= KNOB(KNOB_MAX_INST_COUNT)->getValue())) ||  (sim_end_condition))) { 
        // please complete sim_end_condition 
        print_heartbeat(); 
        print_stats();
        return TRUE; 
    }
    cycle_count++; 
    if (!(cycle_count%5000)) { //print_heartbeat every 5000 cycles
      print_heartbeat(); 
    }


    main_memory->run_a_cycle();        // *NEW-LAB2 
    WB_stage(); 
    MEM_stage(main_memory);  // please modify MEM_stage function argument  /** NEW-LAB2 */
    EX_stage(); 
    ID_stage(); 
    FE_stage(); 
    if (KNOB(KNOB_PRINT_PIPE_FREQ)->getValue() && !(cycle_count%KNOB(KNOB_PRINT_PIPE_FREQ)->getValue())) print_pipeline();
  }
  return TRUE; 
}


/*******************************************************************/
/* Complete the following functions.  */
/* You can add new data structures and also new elements to Op, Pipeline_latch data structure */ 
/*******************************************************************/

void init_structures(memory_c *main_memory) // please modify init_structures function argument  /** NEW-LAB2 */ 
{
  init_op_pool(); 
  init_op_latency();
  init_latches();
  init_registers(); // initializing processor registers to TRUE since no data dependency when processor is powered on
  main_memory->init_mem(); // initializing DRAM and MSHR
  cache_init(data_cache, KNOB(KNOB_DCACHE_SIZE)->getValue(), KNOB(KNOB_BLOCK_SIZE)->getValue(), KNOB(KNOB_DCACHE_WAY)->getValue(), "L1 Cache"); // initializing cache
  use_bpred= KNOB(KNOB_USE_BPRED)->getValue();
  use_vmem= KNOB(KNOB_ENABLE_VMEM)->getValue();
  branchpred = bpred_new(BPRED_GSHARE, (KNOB(KNOB_BPRED_HIST_LEN)->getValue())); // initializing branch_predictor
  dtlb = tlb_new((KNOB(KNOB_TLB_ENTRIES)->getValue()));        // initializing TLB
  ghr_len= (KNOB(KNOB_BPRED_HIST_LEN)->getValue());
  dummy_op=new Op();
  init_op(dummy_op);
}




void WB_stage()
{   
	Op *wb_op;

	if(MEM_latch->op != NULL)
	{
	  if(finalEnd && EX_latch->op->last_op==true && last_op_in_queue==true) //sim_end_condition
	  {
	    	  sim_end_condition=TRUE;
	       
	  }
	}
	
    if(MEM_latch->op_valid) 
    { 
   
     wb_op=MEM_latch->op; // copying MEM_latch to op structure

     if(!use_bpred)
     {
    	 if(wb_op->cf_type >= CF_BR) // Branch instrn at WB stage removes control dependency from all preceeding stages
    		 cd_stall=FALSE;
     }
     else
     {
    	 if(wb_op->cf_type == CF_CBR) // Conditional Branch instrn at WB stage removes control dependency from FE stage
    	 {
    		 if(cd_stall_addr == wb_op->inst_id)  // If this instrn put the cd_stall, only then remove it
    		   cd_stall=FALSE;

    	 }

     }
    // check for dummy op coming back from mshr
 
    if((int)wb_op->dst==-1) // instruction doesn't write to register
     { 
      
      if(!wb_op->last_op)
         retired_instruction++;
      MEM_latch->op_valid=FALSE;
      free_op(wb_op);
     }	
    else
     { 
      register_file[(int)wb_op->dst].valid=TRUE; // to remove register(data) dependency at ID stage

      // Check if the write was an integer write or FP write for power simulations

      if((wb_op->opcode >=2 && wb_op->opcode <=4) || (wb_op->opcode >=7 && wb_op->opcode <=14))
    	  int_regfile_writes++;

      dd_stall=FALSE;  
      if(!wb_op->last_op)
           retired_instruction++;
      MEM_latch->op_valid=FALSE;
      free_op(wb_op);
     }
    }
   
    // Incrementing count for FP instruction reg write for power simulation separately
    // because of malformed FP ops in trace file
    if(wb_op->opcode >=15 && wb_op->opcode <=23)
       float_regfile_writes++;
  
  // Also retire instructions from the queue coming from DRAM queue
    
  while(!mem_op_queue.empty())
    { 
	  wb_op = mem_op_queue.front();

	  if(!use_bpred)
	  {
	   if(wb_op->cf_type >= CF_BR) // Branch instrn at WB stage removes control dependency from all preceeding stages
	  	 cd_stall=FALSE;
	  }
	  else
	  {
	   if(wb_op->cf_type == CF_CBR) // Conditional Branch instrn at WB stage removes control dependency from FE stage
	  	  {
		   if(cd_stall_addr==wb_op->inst_id)
		      cd_stall=FALSE;

	  	  }
	  }

	  
	  if((int)wb_op->dst==-1) // instruction doesn't write to register
	  {
		if(!wb_op->last_op)
		  retired_instruction++;
	    free_op(wb_op);
	  }
	  else
	  {
	    register_file[(int)wb_op->dst].valid=TRUE; // to remove register(data) dependency at ID stage

	    // Check if the write was an integer write or FP write for power simulations

	    if((wb_op->opcode >=2 && wb_op->opcode <=4) || (wb_op->opcode >=7 && wb_op->opcode <=14))
	       	  int_regfile_writes++;

	    dd_stall=FALSE;
	    if(!wb_op->last_op)
	      retired_instruction++;
	    free_op(wb_op);
	  }
	  
	  // Incrementing count for FP instruction reg write separately because of malformed FP ops in trace file
	  if(wb_op->opcode >=15 && wb_op->opcode <=23)
	      float_regfile_writes++;

	  mem_op_queue.pop_front();
	  
    }  
}




void MEM_stage(memory_c *main_memory)
{  
   bool cache_hit, can_forward_sl, can_forward_ss, piggyback_done, inserted_in_MSHR;
   ADDRINT addr;
   Op* op;
   while(!op_queue.empty())
   {
	   op=op_queue.front();
	   mem_op_queue.push_back(op);
	   op_queue.pop_front();
   }
   
   // When last_op passes from EX_latch to MEM_latch, finalEnd becomes TRUE. Now we should just make sure
   // that MSHR becomes empty to make last_op_in_queue to TRUE and ultimately set the sim_end_condition
   if(finalEnd==true)
   {   
	   if(main_memory->m_mshr.empty())
	     last_op_in_queue=true;
	     
   }
 
   if(!use_vmem) // retain code for lab2
   {
	   if(EX_latch->op_valid) // if instrn in EX_latch can be executed in current cycle
	   {
	     Op *mem_op;


	     mem_op=EX_latch->op; // copying EX_latch to mem_op structure

	     if(mem_op->last_op==TRUE)
	     {
	    	 finalEnd=true;

	     	 // Traverse thru the entire MSHR. If its empty, then set last_op_in_queue to true
	    	 if(main_memory->m_mshr.empty())
	    	    last_op_in_queue=true;

	     }


	     if(mem_op->mem_type > 0) // memory instruction
	     {
	    	 /* Accessing dcache and TLB in slightly different manner
	    	  * Waiting for dcache access latency to become 0 and then
	    	  * accessing TLB to get the VPN to PFN mapping. Actual order
	    	  * should be reversed but for timing simulation, order doesn't
	    	  * matter so much
	    	  */
	    	 if(mem_latency_flag==FALSE && mshr_stall==FALSE)
	    	 {
	    		 dcache_latency=KNOB(KNOB_DCACHE_LATENCY)->getValue();
	    		 mem_latency_flag=TRUE;
	    		 mem_latency_stall=TRUE;

	    	 }

	    	 if(mshr_stall==FALSE)
	    	    dcache_latency--;
	    	 if(dcache_latency==0)
	    	 {
	    		 mem_latency_stall=FALSE;
	    		 mem_latency_flag=FALSE;


	    	 // Check for cache hit or miss only when dcache_latency = 0

	    	 if(mem_op->mem_type==MEM_LD)
	    		 addr=mem_op->ld_vaddr;

	    	 if(mem_op->mem_type==MEM_ST)
	    	     addr=mem_op->st_vaddr;



	    	 if(dcache_access(addr))
	    	 {
	    		 // dcache hit for data. Increment LD and ST counters for power evaluation
	    		 if(mem_op->mem_type==MEM_LD)
	    		    dcache_hit_load++;
	    		 if(mem_op->mem_type==MEM_ST)
	    			dcache_hit_store++;

	    		 dcache_hit_count++;
	    		 MEM_latch->op=mem_op; // move instruction to MEM_latch
	    		 MEM_latch->op_valid=TRUE;
	    		 EX_latch->op_valid=FALSE;
	    		 return;
	    	 }
	    	 else  // cache miss
	    	 {
	    		 /*Since its a cache miss, check if load-store or store-store forwarding is possible. If not, then check if piggybacking is possible.
	    		 If not then insert the mem req in MSHR. If no space in MSHR then stall the pipeline. */

	    		 dcache_miss_count++;

	    		 // check if store-load forwarding is possible
	    		 if(mem_op->mem_type==MEM_LD)
	    		 {
	    			 can_forward_sl=main_memory->store_load_forwarding(mem_op);
	    			 if(can_forward_sl)
	    			 {
	    			 	 store_load_forwarding_count++;
	    			     MEM_latch->op=mem_op; // move instruction to MEM_latch
	    			     MEM_latch->op_valid=TRUE;
	    			     EX_latch->op_valid=FALSE;
	    			     return;
	    			 }

	    		 }

	    		 // check if store-store forarding is possible
	    		 if(mem_op->mem_type==MEM_ST)
	    		 {
	    		 	 can_forward_ss=main_memory->store_store_forwarding(mem_op);
	    		     if(can_forward_ss)
	    		     {

	    		         MEM_latch->op=mem_op; // move instruction to MEM_latch
	    		         MEM_latch->op_valid=TRUE;
	    		         EX_latch->op_valid=FALSE;
	    		         return;
	    		     }

	    		 }

	    		// check if piggybacking is possible
	    		piggyback_done=main_memory->check_piggyback(mem_op);
	    		if(piggyback_done)
	    		{
	    		  EX_latch->op_valid=FALSE;
	    		  return;
	    		}
	    		else // mem req can't be piggybacked
	    		{
	    			//insert into MSHR
	    			inserted_in_MSHR=main_memory->insert_mshr(mem_op);
	    			if(!inserted_in_MSHR) // MSHR is full
	    			{
	    				mshr_stall=TRUE;
	    				return;
	    			}
	    			else
	    			{
	    				EX_latch->op_valid=FALSE;
	    				mshr_stall=FALSE;
	    				last_op_in_queue=false;
	    				return;
	    			}
	    		}
	    	 }
	     }
	     }

	     else // not a memory instruction, so move it to the MEM_latch directly
	     {
	     MEM_latch->op=mem_op; // move instruction to MEM_latch
	     MEM_latch->op_valid=TRUE;
	     EX_latch->op_valid=FALSE;
	     return;
		 }
	   }
   }
   else // use virtual memory
   {
	   if(EX_latch->op_valid) // if instrn in EX_latch can be executed in current cycle
	   {
	     Op *mem_op;


	     mem_op=EX_latch->op; // copying EX_latch to mem_op structure

	     if(mem_op->last_op==TRUE)
	     {
	    	 finalEnd=true;

	     	 // Traverse thru the entire MSHR. If its empty, then set last_op_in_queue to true
	    	 if(main_memory->m_mshr.empty())
	    	    last_op_in_queue=true;

	     }


	     if(mem_op->mem_type > 0) // memory instruction
	     {

	    	 vpn_page_size=(KNOB(KNOB_VMEM_PAGE_SIZE)->getValue());
	    	 int threadid=0;
	    	 bool in_tlb;

	    	 // mshr1=2_stall is set when mshr is full while inserting valid mem_op

	    	 if(mshr2_stall)
	    	 {
	    		 inserted_in_MSHR=main_memory->insert_mshr(mem_op);
	    		 if(!inserted_in_MSHR) // MSHR is full
	    		 {
	    		 	mshr2_stall=TRUE;
	    		 	return;
	    		 }
	    		 else
	    		 {
	    		    EX_latch->op_valid=FALSE;
	    		 	mshr2_stall=FALSE;
	    		 	last_op_in_queue=false;
	    		 	return;
	    		 }
	    	 }



	    	 // mshr1_stall is set when mshr is full while inserting dummy_op

	    	 if(mshr1_stall)
	    	 {
	    		 inserted_in_MSHR=main_memory->insert_mshr(dummy_op);
	    		 if(!inserted_in_MSHR) // MSHR is full
	    		 {  mshr1_stall=TRUE;
	    		  	return;
	    		 }
	    		 else
	    		 {	 tlb_stall=TRUE; // stall previous stages because TLB is blocking
	    		 	 mshr1_stall=FALSE;
	    		 	 return;
	    		 }
	    	 }



	    	 if(mem_latency_flag==FALSE && mshr_stall==FALSE && tlb_stall==FALSE)
	    	 {
	    		 dcache_latency=KNOB(KNOB_DCACHE_LATENCY)->getValue();
	    		 mem_latency_flag=TRUE;
	    		 mem_latency_stall=TRUE;

	    	 }

	    	 if(mshr_stall==FALSE && mem_latency_flag==TRUE)
	    	    dcache_latency--;

	    	 if(tlb_stall==TRUE)
	    	 {
	    		 // memory request has come back from DRAM. set dummy_op_flag in fill_queue.
	    		 // Dont send dummy op in broadcast_rdy_op()
	    		 if(dummy_op_flag)
	    		 {
	    			 tlb_stall=FALSE;
	    			 dummy_op_flag=false;
	    			 pfn=vmem_vpn_to_pfn(vpn, threadid);
	    			 tlb_install(dtlb,vpn,threadid, pfn);
	    			 tlb_inst=true; // set to true whenever installed mapping in tlb
	    			 return; // so that TLB entry which is installed will be accessed in the next cycle.
	    		 }
	    	 }

	    	 if(dcache_latency==0 && tlb_stall==FALSE)
	    	 {
	    		 mem_latency_stall=FALSE;
	    		 mem_latency_flag=FALSE;


	    	 // Check for cache hit or miss only when dcache_latency = 0

	    	 if(mem_op->mem_type==MEM_LD)
	    		addr=mem_op->ld_vaddr;

	    	 if(mem_op->mem_type==MEM_ST)
	    	    addr=mem_op->st_vaddr;

	    	 vpn=addr/vpn_page_size;
	    	 in_tlb=tlb_access(dtlb,vpn,threadid, &pfn);

	    	 if(in_tlb)
	    	 {
	    		// TLB hit
	    	   if(tlb_inst==false)
	    	   {
	    		   dtlb_hit_count++;
	    	   }
	    	   paddr=(pfn * vpn_page_size)+(addr % vpn_page_size); // physical address to be used later
	    	   tlb_inst=false;
	    	 }
	    	 else // TLB miss
	    	 {
	    	     dtlb_miss_count++;
	    		 pte=vmem_get_pteaddr(vpn,threadid);
	    		 if(!dcache_access(pte))
	    		 {
	    			 // create a dummy op

	    			 dcache_miss_count++;
	    			 dummy_op->ld_vaddr=pte;
	    			 dummy_op->mem_type=MEM_LD;
	    			 dummy_op->dummy=true;



	    			 //insert into MSHR
	    			 inserted_in_MSHR=main_memory->insert_mshr(dummy_op);
	    			 if(!inserted_in_MSHR) // MSHR is full
	    			 {
	    			 	mshr1_stall=TRUE;
	    			 	return;
	    			 }
	    			 else
	    			 {
	    				 tlb_stall=TRUE; // stall previous stages because TLB is blocking
	    				 mshr1_stall=FALSE;
	    				 return;
	    			 }
	    		 }
	    		 else // PTE entry found in dcache
	    		 {
	    			 // dcache hit for PTE. Increment LD and ST counters for power evaluation
	    			 if(mem_op->mem_type==MEM_LD)
	    			 	    dcache_hit_load++;
	    			 if(mem_op->mem_type==MEM_ST)
	    			 		dcache_hit_store++;

	    			 // Get PFN from vmem_vpn_to_pfn() and call tlb_install()
	    			 // Use this PFN to calculate paddr and continue with dcache access lab2 code below
	    			 dcache_hit_count++;
	    			 pfn=vmem_vpn_to_pfn(vpn, threadid);
	    			 tlb_install(dtlb,vpn,threadid, pfn);
	    			 tlb_inst=true; // set to true whenever installed mapping in tlb
	    			 return;

	    		 }

	    	 }

	    	 // Old code from lab 2

	    	 if(dcache_access(paddr))
	    	 {
	    		 // dcache hit for data. Increment LD and ST counters for power evaluation
	    		 if(mem_op->mem_type==MEM_LD)
	    		     dcache_hit_load++;
	    		 if(mem_op->mem_type==MEM_ST)
	    		 	dcache_hit_store++;

	    		 dcache_hit_count++;
	    		 MEM_latch->op=mem_op; // move instruction to MEM_latch
	    		 MEM_latch->op_valid=TRUE;
	    		 EX_latch->op_valid=FALSE;
	    		 return;
	    	 }
	    	 else  // cache miss
	    	 {
	    		 /*Since its a cache miss, check if load-store or store-store forwarding is possible. If not, then check if piggybacking is possible.
	    		 If not then insert the mem req in MSHR. If no space in MSHR then stall the pipeline. */

	    		 dcache_miss_count++;

	    		 // check if store-load forwarding is possible
	    		 if(mem_op->mem_type==MEM_LD)
	    		 {
	    			 mem_op->ld_vaddr=paddr;
	    			 can_forward_sl=main_memory->store_load_forwarding(mem_op);
	    			 if(can_forward_sl)
	    			 {
	    			 	 store_load_forwarding_count++;
	    			     MEM_latch->op=mem_op; // move instruction to MEM_latch
	    			     MEM_latch->op_valid=TRUE;
	    			     EX_latch->op_valid=FALSE;
	    			     return;
	    			 }

	    		 }

	    		 // check if store-store forarding is possible
	    		 if(mem_op->mem_type==MEM_ST)
	    		 {
	    			 mem_op->st_vaddr=paddr;
	    			 can_forward_ss=main_memory->store_store_forwarding(mem_op);
	    		     if(can_forward_ss)
	    		     {

	    		         MEM_latch->op=mem_op; // move instruction to MEM_latch
	    		         MEM_latch->op_valid=TRUE;
	    		         EX_latch->op_valid=FALSE;
	    		         return;
	    		     }

	    		 }

	    		// check if piggybacking is possible

	    		piggyback_done=main_memory->check_piggyback(mem_op);
	    		if(piggyback_done)
	    		{
	    		  EX_latch->op_valid=FALSE;
	    		  return;
	    		}
	    		else // mem req can't be piggybacked
	    		{
	    			//insert into MSHR
	    			inserted_in_MSHR=main_memory->insert_mshr(mem_op);
	    			if(!inserted_in_MSHR) // MSHR is full
	    			{
	    				mshr2_stall=TRUE;
	    				return;
	    			}
	    			else
	    			{
	    				EX_latch->op_valid=FALSE;
	    				mshr2_stall=FALSE;
	    				last_op_in_queue=false;
	    				return;
	    			}
	    		}
	    	 }
	     }
	     }

	     else // not a memory instruction, so move it to the MEM_latch directly
	     {
	     MEM_latch->op=mem_op; // move instruction to MEM_latch
	     MEM_latch->op_valid=TRUE;
	     EX_latch->op_valid=FALSE;
	     return;
		 }
	   }
   }

}





void EX_stage()
{
   /* The very first execution(for any instrn) in EX_stage would always match the 2nd condition.
      If latency>1 for an instrn, then 1st condition is matched in subsequent cycles in EX stage(so decrement latency)
      Instrn moves to EX latch when latency becomes 0 */
   
  // If any memory stage stall=TRUE, make sure EX_latch valid bit is TRUE so that MEM_stage executes in next cycle. 

   if(ex_stall==TRUE)
     {
      Op *ex_op=ID_latch->op; // read from ID_latch
      latency--;
      if(latency==0)
       {
    	  ex_stall=FALSE;
    	  EX_latch->op=ex_op; 
    	  EX_latch->op_valid=TRUE;
    	  ID_latch->op_valid=FALSE;
    	  if(ex_op->last_op==true)
    		  finalEnd=true;
       }
       return;
     }

   if(use_vmem || use_bpred)
   {
	   if(mem_latency_stall==FALSE && mshr_stall==FALSE && mshr1_stall==FALSE && mshr2_stall==FALSE && ex_stall==FALSE && tlb_stall==FALSE && tlb_inst==false) // lab3 specific stalls
	         {

	           if(ID_latch->op_valid==TRUE)
	           {
	            Op *ex_op=ID_latch->op; // read from ID_latch
	            latency=get_op_latency(ex_op);

	            if(latency>1)
	            	{
	         	   ex_stall=TRUE;
	         	   latency--;
	         	   return;
	            	}
	            else // single cycle instruction, so pass it to EX latch
	            {
	         	   EX_latch->op=ex_op;
	         	   EX_latch->op_valid=TRUE;
	         	   ID_latch->op_valid=FALSE;
	         	   if(ex_op->last_op==true)
	         	       	finalEnd=true;
	         	   return;
	            }
	           }
	         }
   }
   else // lab2 code retained
   {
	   if(mem_latency_stall==FALSE && mshr_stall==FALSE && ex_stall==FALSE)
	         {

	           if(ID_latch->op_valid==TRUE)
	           {
	            Op *ex_op=ID_latch->op; // read from ID_latch
	            latency=get_op_latency(ex_op);

	            if(latency>1)
	            	{
	         	   ex_stall=TRUE;
	         	   latency--;
	         	   return;
	            	}
	            else // single cycle instruction, so pass it to EX latch
	            {
	         	   EX_latch->op=ex_op;
	         	   EX_latch->op_valid=TRUE;
	         	   ID_latch->op_valid=FALSE;
	         	   if(ex_op->last_op==true)
	         	       	finalEnd=true;
	         	   return;
	            }
	           }
	         }
   }
    

}  



  
void ID_stage()
{  
  /* If there is an ex_stall or mem_latency_stall in the current cycle, ID stage should not do any work. 
     If there is dd_stall, increment data_hazard_count but don't do any other work.
     If none of the stalls are in place, then enter ID stage only if FE latch's valid flag is TRUE. */

	if(!use_bpred) // retain lab2 code
	{
		if(dd_stall==TRUE)
		  data_hazard_count++;

		if(FE_latch->op_valid==TRUE)
		{
			Op *id_op;
			id_op=FE_latch->op;
			// Count number of scheduler accesses for Integer and FP instructions for power simulation

			if(id_op->opcode >=15 && id_op->opcode <=23)
			   fp_sched_access++;
			if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
		       int_sched_access++;
		}


		  if(mem_latency_stall==FALSE && mshr_stall==FALSE && mshr1_stall==FALSE && mshr2_stall==FALSE && dd_stall==FALSE && ex_stall==FALSE && tlb_stall==FALSE && tlb_inst==false)
		  {
		   if(FE_latch->op_valid==TRUE)
		   {

		    Op *id_op;
		    id_op=FE_latch->op;

		   // check for control dependency

		   if(id_op->cf_type >= CF_BR)
		      cd=TRUE;
		   else
		      cd=FALSE;


		  // check for data dependency
		  int src0_reg=(int)id_op->src[0];
		  int src1_reg=(int)id_op->src[1];

		  int dst_reg=(int)id_op->dst;



		  if(src0_reg==-1 && src1_reg!=-1)  //instruction has src1 register only
		   {
		     if(register_file[src1_reg].valid==FALSE) // src1 register is not ready
		     dd=TRUE;
		     else
		     dd=FALSE;

		     // Count number of register reads for power calculation

		    if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
		         int_reg_read++;
		   }



		   if(src1_reg==-1 && src0_reg!=-1)  //instruction has src0 register only
		   {
		     if(register_file[src0_reg].valid==FALSE) // src0 register is not ready
		      dd=TRUE;
		     else
		      dd=FALSE;

		     // Count number of register reads for power calculation

		    if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
		           int_reg_read++;

		   }


		   if(src1_reg!=-1 && src0_reg!=-1) ////instruction has both src0 and src1 registers
		   {
		     if(register_file[src0_reg].valid==FALSE || register_file[src1_reg].valid==FALSE) // either register is not ready
		     dd=TRUE;
		     else
		     dd=FALSE;

		     // Count number of register reads for power calculation

		    if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
		           int_reg_read+=2;
		   }

		   // Assuming every FP instruction to have 2 src registers because of malformed FP ops in trace file
		    if(id_op->opcode >=15 && id_op->opcode <=23)
			      	   fp_reg_read+=2;

		   if(src1_reg==-1 && src0_reg==-1)  // if there are no src registers, then no data dependency
		      dd=FALSE;


		    /* if CD and DD => stall FE for current cycle. Dont pass current instruction to ID latch. CH++ for every cycle. DH++
		    if CD but no DD=> stall FE for current cycle. Pass current instruction to ID latch. CH++ for just one cycle.
		    if no CD but DD=> stall FE. DH++. dont pass current instruction to ID latch
		    if no CD and no DD=> no stalls. pass the current instruction to ID latch     */

		    if(cd==TRUE && dd==TRUE)
		     {
		       cd_stall=TRUE;
		       dd_stall=TRUE;
		       control_hazard_count++;
		       data_hazard_count++;
		       return;
		     }

		    if(cd==TRUE && dd==FALSE)
		     {
		       dd_stall=FALSE;
		       cd_stall=TRUE;
		       control_hazard_count++;
		       ID_latch->op=id_op;
		       ID_latch->op_valid=TRUE;
		       FE_latch->op_valid=FALSE;
		       if(dst_reg!=-1) // instruction writes to destination register, so change valid flag of dst to FALSE
		         register_file[dst_reg].valid=FALSE;
		       return;
		     }

		    if(cd==FALSE && dd==TRUE)
		     {
		       cd_stall=FALSE;
		       dd_stall=TRUE;
		       data_hazard_count++;
		       return;
		     }

		    if(cd==FALSE && dd==FALSE)
		     {
		       cd_stall=FALSE;
		       dd_stall=FALSE;
		       ID_latch->op=id_op;
		       ID_latch->op_valid=TRUE;
		       FE_latch->op_valid=FALSE;
		       if(dst_reg!=-1) // instruction writes to destination register, so change valid flag to FALSE
		         register_file[dst_reg].valid=FALSE;
		       return;
		     }

		  }

		     }
	}
	
	else // bpred is used
	{
		if(dd_stall==TRUE)
		  data_hazard_count++;

		if(FE_latch->op_valid==TRUE)
		{
			Op *id_op;
			id_op=FE_latch->op;
			// Count number of scheduler accesses for Integer and FP instructions

			if(id_op->opcode >=15 && id_op->opcode <=23)
			   fp_sched_access++;
			if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
			   int_sched_access++;
		}

		if(mem_latency_stall==FALSE && mshr_stall==FALSE && mshr1_stall==FALSE && mshr2_stall==FALSE && dd_stall==FALSE && ex_stall==FALSE && tlb_stall==FALSE && tlb_inst==false)
		{

		 if(FE_latch->op_valid==TRUE)
		 {

		   Op *id_op;
		   id_op=FE_latch->op;

		  // control dependency check is not required when BPRED is enabled


		  // check for data dependency
		  int src0_reg=(int)id_op->src[0];
		  int src1_reg=(int)id_op->src[1];

		  int dst_reg=(int)id_op->dst;



		  if(src0_reg==-1 && src1_reg!=-1)  //instruction has src1 register only
		   {
		     if(register_file[src1_reg].valid==FALSE) // src1 register is not ready
		     dd=TRUE;
		     else
		     dd=FALSE;

		     // Count number of register reads for power calculation

		    if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
		    	int_reg_read++;
		   }



		   if(src1_reg==-1 && src0_reg!=-1)  //instruction has src0 register only
		   {
		     if(register_file[src0_reg].valid==FALSE) // src0 register is not ready
		      dd=TRUE;
		     else
		      dd=FALSE;

		     // Count number of register reads for power calculation

		    if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
		          int_reg_read++;
		   }


		   if(src1_reg!=-1 && src0_reg!=-1) ////instruction has both src0 and src1 registers
		   {
		     if(register_file[src0_reg].valid==FALSE || register_file[src1_reg].valid==FALSE) // either register is not ready
		     dd=TRUE;
		     else
		     dd=FALSE;

		     // Count number of register reads for power calculation

		    if((id_op->opcode >=2 && id_op->opcode <=4) || (id_op->opcode >=7 && id_op->opcode <=14))
		          int_reg_read++;
		   }

		   // Assuming every FP instruction to have 2 src registers because of malformed FP ops in trace file
		   if(id_op->opcode >=15 && id_op->opcode <=23)
		       	   fp_reg_read+=2;

		   if(src1_reg==-1 && src0_reg==-1)  // if there are no src registers, then no data dependency
		      dd=FALSE;


		    /*
		    if no CD but DD=> stall FE. DH++. dont pass current instruction to ID latch
		    if no CD and no DD=> no stalls. pass the current instruction to ID latch     */

		    if(dd==TRUE)
		     {
		       dd_stall=TRUE;
		       data_hazard_count++;
		       return;
		     }

		    if(dd==FALSE)
		     {
		       dd_stall=FALSE;
		       ID_latch->op=id_op;
		       ID_latch->op_valid=TRUE;
		       FE_latch->op_valid=FALSE;
		       if(dst_reg!=-1) // instruction writes to destination register, so change valid flag to FALSE
		         register_file[dst_reg].valid=FALSE;
		       return;
		     }

		  }

		}
	}
}
 



void FE_stage()
{
  /* FE latch valid flag has to be FALSE only when there is a cd_stall AND not other stalls. 
     Because only in a cd_stall, ID stage has to continue processing the instructions and move them
     to the next latch. Whereas when there are other stalls(ex_stall, dd_stall etc.) alongwith cd_stall, 
     ID stage is also stalled. */
	
	 fe_ctr++;
   if(fe_ctr>100000) // If simulator waits for an op for > 1000 cycles, terminate the program
    {
     //std::cout<<"****************Detected a loop. Terminating the simulator***************"<<endl;
     exit(1);
    } 
   

   if(tlb_stall==FALSE && mem_latency_stall==FALSE && mshr_stall==FALSE && mshr1_stall==FALSE && mshr2_stall==FALSE && cd_stall==FALSE && dd_stall==FALSE && ex_stall==FALSE && fe_stall==FALSE && tlb_inst==false)  // added stalls in lab2 code
   {

   
     if(FE_latch->op_valid==FALSE) // don't overwrite instrn in FE latch since ID stage has not used it yet
     {
  
 
      Op *fe_op = get_free_op(); // op currently in FE_stage
      bool success=get_op(fe_op);

      // Count different types of instructions
      if(fe_op->opcode >=15 && fe_op->opcode <=23)
    	  fp_instruction++;
      if((fe_op->opcode >=2 && fe_op->opcode <=4) || (fe_op->opcode >=7 && fe_op->opcode <=14))
    	  int_instruction++;

      if(fe_op->opcode == OP_CF)
    	  branch_instruction++;
      if(fe_op->opcode == OP_LD)
    	  load_instruction++;
      if(fe_op->opcode == OP_ST)
    	  store_instruction++;

      if(fe_op->opcode == OP_IMUL || fe_op->opcode == OP_MM)
    	  mul_instruction++;



      if(use_bpred) // perform branch prediction only if asked through CLI
      {
    	  if(fe_op->cf_type==CF_CBR) // only for conditional branch instrn, predict branch direction
    	  {
    	   	btaken=bpred_access(branchpred, fe_op->instruction_addr); // btaken is the predicted branch direction
    	   	if(btaken != fe_op->actually_taken)
    	   	{
    	   		bpred_mispred_count++;
    	   		cd_stall= TRUE; // mispredicted branch, so stall the FE stage for subsequent cycles till branch is resolved
    	   		cd_stall_addr=fe_op->inst_id;
    	   	}
    	   	else
    	   		bpred_okpred_count++;

    	   	bpred_update(branchpred, fe_op->instruction_addr, (int)btaken, (int)fe_op->actually_taken);
          }

      }

   
      if(success==FALSE) //last instruction
      {  
    	  fe_op->last_op=TRUE;
    	  fe_stall=TRUE;

      }
      
      
  
      FE_latch->op=fe_op;
      FE_latch->op_valid=TRUE; 
      fe_ctr=0; // Instruction was successfully read, so reset the loop detection counter
     }
   } 
  //   next_pc = pc + op->inst_size;  // you need this code for building a branch predictor 
   if(cd_stall)
      {
   	   control_hazard_count++;
      }
}


void  init_latches()
{
  MEM_latch = new pipeline_latch();
  EX_latch = new pipeline_latch();
  ID_latch = new pipeline_latch();
  FE_latch = new pipeline_latch();

  MEM_latch->op = NULL;
  EX_latch->op = NULL;
  ID_latch->op = NULL;
  FE_latch->op = NULL;

  /* you must set valid value correctly  */ 
  MEM_latch->op_valid = false;
  EX_latch->op_valid = false;
  ID_latch->op_valid = false;
  FE_latch->op_valid = false;

}


bool icache_access(ADDRINT addr) 
{   /** please change uint32_t to ADDRINT NEW-LAB2 */ 
  /* For Lab #1, you assume that all I-cache hit */     
  bool hit = FALSE; 
  if (KNOB(KNOB_PERFECT_ICACHE)->getValue()) hit = TRUE; 
  return hit; 
}



bool dcache_access(ADDRINT addr) 
{ /** please change uint32_t to ADDRINT NEW-LAB2 */ 
  /* For Lab #1, you assume that all D-cache hit */     
  /* For Lab #2, you need to connect cache here */   // NEW-LAB2  

  bool hit = FALSE;
  if (KNOB(KNOB_PERFECT_DCACHE)->getValue()) 
	  hit = TRUE; 
  else
	  hit=(bool)cache_access(data_cache,addr);

  return hit; 
}



// NEW-LAB2 
void dcache_insert(ADDRINT addr)  // NEW-LAB2 
{                                 
  /* dcache insert function */   
  cache_insert(data_cache, addr) ;
  inserted_in_dcache++; // Counter for power evaluation
 
}                                       



void broadcast_rdy_op(Op* op)              
{                                          
  // mem ops are done.  move the op into WB stage   // NEW-LAB2 
  // Queue the op into a queue to retire multiple ops in one cycle
	if(use_vmem)
	{
		if(op->mem_read_size==0 && op->mem_write_size==0 && op->dummy==true)
		{
			//dummy instruction has come back from DRAM. Don't push it to op_queue for retirement
			dummy_op_flag=true;
			return;
		}
	  op_queue.push_back(op);
	  return;
	}
	
	op_queue.push_back(op);

}      



/* utility functions that you might want to implement */     // NEW-LAB2 
int get_dram_row_id(ADDRINT addr)    // NEW-LAB2 
{  
  // addr >> 6;   // NEW-LAB2 
	
	int dram_page_size,dram_row_id;
	
	dram_page_size = KNOB(KNOB_DRAM_PAGE_SIZE)->getValue();
	dram_row_id=(addr/(dram_page_size*1024));
	
	return dram_row_id;

 // return 2;   
}  



int get_dram_bank_id(ADDRINT addr)  // NEW-LAB2 
{  // NEW-LAB2 
  // (addr >> 6);   // NEW-LAB2 

	int dram_page_size,dram_bank_num,dram_bank_id;
	
	dram_page_size = KNOB(KNOB_DRAM_PAGE_SIZE)->getValue();
	dram_bank_num = KNOB(KNOB_DRAM_BANK_NUM)->getValue();
	dram_bank_id=(addr/(dram_page_size*1024))%dram_bank_num;
	
	return dram_bank_id;
	
  //return 1;   // NEW-LAB2 
}  // NEW-LAB2 


