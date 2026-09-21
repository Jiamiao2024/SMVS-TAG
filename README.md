
# SMVS-TAG

MATLAB code for the paper "Scalable Multi-View Subspace Clustering with Tensorized Anchor Guidance"

## Abstract

Anchor-based multi-view clustering methods have gained significant attention for their effectiveness in handling large-scale datasets in recent years. The performance of these methods is highly dependent on anchor quality.
However, current methods neglect the interactive relationships among cross-view anchors, failing to effectively discover and exploit  consistent and complementary information, leading to noisy or suboptimal anchor representations. In this paper, we propose a novel scalable multi-view subspace clustering method with tensorized anchor guidance, which directly couples anchors across views to improve clustering performance. Specifically, we construct a third-order anchor tensor from view-specific anchors in a low-dimensional latent space. By imposing a tensor Schatten $p$-norm constraint on the anchor tensor, we can explicitly capture cross-view low-rank structure and jointly exploit consistent and complementary information  among anchors. Moreover, the tensorized anchor regularizer is independent of the number of samples, which reduces both time and space complexity.  Experimental results on seven datasets  demonstrate that SMVS-TAG achieves superior effectiveness and stability compared to state-of-the-art large-scale MVC methods.


# Citation

If you find our code useful, please cite:

@inproceedings{jia2026scalable,
  title={Scalable Multi-View Subspace Clustering with Tensorized Anchor Guidance},
  author={Jia, Miao and Hu, Xingchen and Liu, Jiyuan and Wang, Siwei and Wang, Min and Chen, Zijian},
  booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition},
  pages={14367--14376},
  year={2026}
}

Thanks. Any problem can contact Miao Jia (jiamiao0526@163.com).