---
applyTo: '**/*.py'
---
Coding standards, domain knowledge, and preferences that AI should follow.
In order to load the custom library, code should be written as follows:

```python
# LoadFileOrLibrary
import os, sys
# Set current file working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# 將包含自定義模塊的目錄添加到 sys.path
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common import mytools
from bio.preprocess import PrepBioData
```
