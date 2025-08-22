setwd("F:\\新图课题")
rm(list = ls())
library(data.table)
library(tidyr)
library(dplyr)
# 读取RDS文件
pathway_import <- readRDS("C:\\Users\\Admin\\Documents\\WeChat Files\\wxid_iur46pqup1b922\\FileStorage\\File\\2025-05\\pathway_import.rds")
pathway_export <- readRDS("C:\\Users\\Admin\\Documents\\WeChat Files\\wxid_iur46pqup1b922\\FileStorage\\File\\2025-05\\pathway_export.rds")
reaction<-fread('C:\\Users\\Admin\\Desktop\\通路富集\\ReactionDatabase.txt')%>%as.data.frame()


##Import
##函数
enrichment_analysis <- function(pathway_import, diff_metabolites) {
  # 初始化结果数据框
  result <- data.frame(
    pathway = character(),
    p_adj = numeric(),
    count = integer(),
    metabolites = character(),
    stringsAsFactors = FALSE
  )
  
  # 总代谢物数（用于超几何检验）
  total_metabolites <- unique(unlist(pathway_import))
  N <- length(total_metabolites)
  
  # 差异代谢物数
  n <- length(diff_metabolites)
  
  # 遍历每个通路
  for (pathway in names(pathway_import)) {
    # 通路中的代谢物
    pathway_metabolites <- pathway_import[[pathway]]
    K <- length(pathway_metabolites)
    
    # 差异代谢物中属于该通路的
    overlap <- intersect(diff_metabolites, pathway_metabolites)
    k <- length(overlap)
    
    if (k > 0) {
      # 超几何检验 p-value
      p_value <- phyper(
        q = k - 1,  # P(X >= k)
        m = K,
        n = N - K,
        k = n,
        lower.tail = FALSE
      )
      
      # 保留前10个代谢物
      if (k > 10) {
        overlap <- overlap[1:10]
      }
      
      # 添加到结果
      result <- rbind(result, data.frame(
        pathway = pathway,
        p_adj = p_value,  # 这里可以进一步做多重检验校正
        count = k,
        metabolites = paste(overlap, collapse = "/"),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 对 p-value 进行多重检验校正（BH方法）
  result$p_adj <- p.adjust(result$p_adj, method = "BH")
  
  # 按 p_adj 排序
  result <- result[order(result$p_adj), ]
  
  return(result)
}
read.csv("C:\\Users\\Admin\\Desktop\\通路富集\\diff_metabolites.csv")->metabolites
diff<-metabolites$suoxie

diff<-gsub("\\[.*?\\]", "", diff)#删除[]和[]中的内容
use_pathway <- enrichment_analysis(pathway_import, diff)
colnames(use_pathway)<-c('Description','p.adjust','Count','geneID')
use_pathway$Category<-reaction$Subsystem_general[match(use_pathway$Description,reaction$Subsystem)]
use_pathway[-which(use_pathway$Category=="Transport"),]->use_pathway
use_pathway<-use_pathway[1:10,]



width <- 0.3#最左侧色块的宽度、色块与圆点、圆点与y轴之间的宽度
# x 轴长度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
# 左侧分类标签数据
rect.data <- group_by(use_pathway, Category) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

use_pathway$Description<-factor(use_pathway$Description,levels = use_pathway$Description[order(use_pathway$Category)])
library(ggplot2)
#if (!require("devtools", quietly = TRUE))
# install.packages("devtools")

#devtools::install_github("dxsbiocc/gground")
library(gground)
#install.packages("ggprism")
library(ggprism)

use_pathway$title<-'Production'
Import<-use_pathway %>%
  ggplot(aes(-log10(p.adjust), y = Description, fill = Category)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8 #width调整矩形宽度
  ) +
  geom_text(
    aes(x = 0.05,label = Description),
    hjust = 0, size = 5
  ) +
  geom_text(
    aes(x = 0.1,label = geneID, colour = Category), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -width, label = Count)
  ) +
  scale_size_continuous(name = 'Count', 
                        range = c(5, 8),
                        labels = scales::number_format(accuracy = 1)) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~title)+
  labs(y = NULL) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 2), 
    expand = expansion(c(0.01, 0))
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(face = "bold",size = 12)#标题字体大小
  )
Import  











##export
##函数
enrichment_analysis <- function(pathway_export, diff_metabolites) {
  # 初始化结果数据框
  result <- data.frame(
    pathway = character(),
    p_adj = numeric(),
    count = integer(),
    metabolites = character(),
    stringsAsFactors = FALSE
  )
  
  # 总代谢物数（用于超几何检验）
  total_metabolites <- unique(unlist(pathway_export))
  N <- length(total_metabolites)
  
  # 差异代谢物数
  n <- length(diff_metabolites)
  
  # 遍历每个通路
  for (pathway in names(pathway_import)) {
    # 通路中的代谢物
    pathway_metabolites <- pathway_export[[pathway]]
    K <- length(pathway_metabolites)
    
    # 差异代谢物中属于该通路的
    overlap <- intersect(diff_metabolites, pathway_metabolites)
    k <- length(overlap)
    
    if (k > 0) {
      # 超几何检验 p-value
      p_value <- phyper(
        q = k - 1,  # P(X >= k)
        m = K,
        n = N - K,
        k = n,
        lower.tail = FALSE
      )
      
      # 保留前10个代谢物
      if (k > 10) {
        overlap <- overlap[1:10]
      }
      
      # 添加到结果
      result <- rbind(result, data.frame(
        pathway = pathway,
        p_adj = p_value,  # 这里可以进一步做多重检验校正
        count = k,
        metabolites = paste(overlap, collapse = "/"),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 对 p-value 进行多重检验校正（BH方法）
  result$p_adj <- p.adjust(result$p_adj, method = "BH")
  
  # 按 p_adj 排序
  result <- result[order(result$p_adj), ]
  
  return(result)
}
read.csv("C:\\Users\\Admin\\Desktop\\通路富集\\diff_metabolites.csv")->metabolites
diff<-metabolites$suoxie

diff<-gsub("\\[.*?\\]", "", diff)#删除[]和[]中的内容
use_pathway <- enrichment_analysis(pathway_export, diff)
colnames(use_pathway)<-c('Description','p.adjust','Count','geneID')
use_pathway$Category<-reaction$Subsystem_general[match(use_pathway$Description,reaction$Subsystem)]
use_pathway[-which(use_pathway$Category=="Transport"),]->use_pathway
use_pathway<-use_pathway[1:10,]



width <- 0.3#最左侧色块的宽度、色块与圆点、圆点与y轴之间的宽度
# x 轴长度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
# 左侧分类标签数据
rect.data <- group_by(use_pathway, Category) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

use_pathway$Description<-factor(use_pathway$Description,levels = use_pathway$Description[order(use_pathway$Category)])
library(ggplot2)
#if (!require("devtools", quietly = TRUE))
# install.packages("devtools")

#devtools::install_github("dxsbiocc/gground")
library(gground)
#install.packages("ggprism")
library(ggprism)

use_pathway$title<-'Consumption'

Export<-use_pathway %>%
  ggplot(aes(-log10(p.adjust), y = Description, fill = Category)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8 #width调整矩形宽度
  ) +
  geom_text(
    aes(x = 0.05,label = Description),
    hjust = 0, size = 5
  ) +
  geom_text(
    aes(x = 0.1,label = geneID, colour = Category), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -width, label = Count)
  ) +
  scale_size_continuous(name = 'Count', 
                        range = c(5, 8),
                        labels = scales::number_format(accuracy = 1)) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~title)+
  labs(y = NULL) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 2), 
    expand = expansion(c(0.01, 0))
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(face = "bold",size = 12)#标题字体大小
  )
Export  



# 安装cowplot包（如果尚未安装）
if (!require("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
library(cowplot)
# 将Import和Export两个图拼在一起
combined_plot <- plot_grid(Import, Export, ncol = 2)
# 显示组合后的图
print(combined_plot)
ggsave("combined_plot.png", combined_plot, width = 30, height = 6) 
ggsave("combined_plot.pdf", combined_plot, width = 30, height = 6) 
