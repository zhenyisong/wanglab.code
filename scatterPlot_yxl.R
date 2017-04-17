library(outliers)
"
raw data is from Yxl's excel file
the file was sent by Yuan and recieved.
the prefix of her sent file is 8-4*
the following data was manually curated
from her Excel file.
"
ctrl.1 <- list(actin = c(19, 20, 15),thoc5 = c(21, 19, 22))
ctrl.2 <- list(actin = c(14, 17, 19),thoc5 = c(12, 17, 16))
ctrl.3 <- list(actin = c(20, 18, 22),thoc5 = c(24, 26, 20))
ctrl.4 <- list(actin = c(14, 10, 12),thoc5 = c(12, 17, 11))
ctrl.5 <- list(actin = c(11, 19, 15),thoc5 = c(12, 15, 17))
as.1   <- list(actin = c( 7,  9, 11),thoc5 = c( 6, 10, 17))
as.2   <- list(actin = c( 2,  8, 11),thoc5 = c( 7,  3, 14))
as.3   <- list(actin = c( 6,  7,  5),thoc5 = c( 7, 10,  2))
as.4   <- list(actin = c(10,  5,  3),thoc5 = c( 1,  8,  2))
as.5   <- list(actin = c( 3,  4,  2),thoc5 = c( 3, 11, 10))

ctrl_1 <- median(ctrl.1$thoc5)/median(ctrl.1$actin)
ctrl_2 <- median(ctrl.2$thoc5)/median(ctrl.2$actin)
ctrl_3 <- median(ctrl.3$thoc5)/median(ctrl.3$actin)
ctrl_4 <- median(ctrl.4$thoc5)/median(ctrl.4$actin)
ctrl_5 <- median(ctrl.5$thoc5)/median(ctrl.5$actin)

as_1   <- median(as.1$thoc5)/median(as.1$actin)
as_2   <- median(as.2$thoc5)/median(as.2$actin)
as_3   <- median(as.3$thoc5)/median(as.3$actin)
as_4   <- median(as.4$thoc5)/median(as.4$actin)
as_5   <- median(as.5$thoc5)/median(as.5$actin)

control.set <- c(ctrl_1,ctrl_2,ctrl_3,ctrl_4,ctrl_5)
patient.set <- c(as_1,as_2,as_3,as_4,as_5)

shapiro.test(control.set)
shapiro.test(patient.set)

grubbs.test(control.set)
grubbs.test(patient.set)

var.test(control.set,patient.set)

t.test(control.set, patient.set, var.equal = FALSE)

norm.data <- c(control.set, patient.set )
data.grp  <- factor(c(rep(1,5),rep(0,5)), levels = c(1, 0), labels = c('Control', 'as'))


boxplot( norm.data ~ data.grp, ylim = c(0,4),
         pars = list( boxwex    = 0.3, 
                      staplewex = 0.5, 
                      outwex    = 0.5),
         outline = F,
         at   = c(1, 1 + 0.5))
         
stripchart( norm.data ~ data.grp, vertical = T, 
            method = "jitter", pch = 16, 
            at     = c(1, 1 + 0.5),
            offset = 0.3,
            col = c("blue", 'red'),
            add = TRUE)
segments(1, 3.5, 1.5, 3.5)
segments(1, 3.2, 1, 3.5)
segments(1 + 0.5, 3.2, 1 + 0.5, 3.5)
text(1.2, 3.7, 'N.S')

