
set.seed(123)
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Whole_dataset.Rda")
dataset<-dataset[!is.na(dataset$SYNTAX.SCORE),]
train_indices <- createMultiFolds(
  y     = dataset$Binary.SS,
  k     = 10,
  times = 10
)
n <- length(dataset$Binary.SS)

test_indices <- lapply(
  train_indices,
  function(train_idx) setdiff(seq_len(n), train_idx)
)

test_indices_all<-test_indices
for(outer in 1:10){
  test_indices<-list()
  for(inner in 1:10){
    name_fold <- sprintf(
      "Fold%02d.Rep%02d",
      inner,
      outer
    )
    
    test_indices[[inner]]<-test_indices_all[[name_fold]]
  }
  name_file<-sprintf("C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold%01d/Test_inds.Rda",outer)
  save(test_indices,file=name_file)
}




