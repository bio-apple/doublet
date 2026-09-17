
1.我在当前文件夹下载一些文献，你参考这些文献，帮我写一个单细胞数据分析的多胞率分析的skill.md， 另外参考来源还有https://bioconductor.org/books/3.23/OSCA.advanced/doublet-detection.html，要求针对同一个输入数据，给出多种算法的预测结果

2. 算法选择
    - scDblFinder
    - Scrublet
    - cxds
    - bcds
    - hybrid
    - Solo
    - DoubletDetection
    - DoubletFinder
    - DoubletDecon

3. 只做 RNA-only

4. 输入大于2万个细胞不使用对时间敏感的算法

5. 永远不把 expected rate 当 cutoff 输入。只保留不依赖该参数就能出的 Score；若要 Call，只用数据驱动阈值（Scrublet 模拟分数双峰、scDblFinder dbr.sd=1、cxds/bcds 的 score 分布/outlier）。nExp 驱动的 DoubletFinder 二元判定仍视为违规

6. 两者都交：Score（必须）+ 数据驱动的 Call + 由 Call 算出的 Predicted Rate；不移除细胞