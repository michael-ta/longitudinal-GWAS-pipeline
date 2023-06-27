
process GALLOPCOX_INPUT {

  scratch true
  label 'small'
  
  input:
    tuple val(chrname), path(plink_input) //from p3_get_plink_chunks_pfile_input
  output:
    tuple val(chrname), path("allchr_${params.dataset}_p2in.txt") //into p3_plink_chunks
    
  // when:
  //  params.longitudinal_flag || params.survival_flag
  
  script:
    """
    #!/usr/bin/env python3
    fn = "${plink_input}"
    out_fn = "allchr_${params.dataset}_p2in.txt"
    count = 0
    id_pairs = []
    start, end =  None,None
    with open(fn, 'r') as f:
      for l in iter(f.readline, ''):
        if l[0] == '#':
          continue
        data = l.strip().split('\t')
        vid = data[2]
        count += 1
        if start is None:
          start = vid

        if count >= ${params.plink_chunk_size}:
          end = vid
          id_pairs.append( (start, end) )
          start = None
          end = None
          count = 0

    if count > 0:
      end = vid
      id_pairs.append( (start, end) )

    with open(out_fn, 'w') as f:
      for start,end in id_pairs:
        f.write( '\\t'.join([start,end]) + '\\n' )
    """
}