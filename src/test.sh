function test_pipeline1 {
  cmd='python -c "import time;print('Pipeline1:');time.sleep(3)"'
  jobname="test_pipeline1"
}

function test_pipeline2 {
  cmd='python -c "import time;print('Pipeline2:');time.sleep(3)"'
  jobname="test_pipeline2"
}

function test_pipeline3 {
  cmd='python -c "import time;print('Pipeline3:');time.sleep(60)"'
  jobname="test_pipeline3"
}
