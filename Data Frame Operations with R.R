library(dplyr)
library(reshape2)
titanic <- read.csv("https://raw.githubusercontent.com/hurratanvir/algo-trade/master/titanic.csv")

#The titanic dataset contains the data about the passengers 
#(Names, Fare etc. and whether they survived or not) 
#who boarded the famous titanic ship in 1912.
#The various columns in the data frame are as:

#PassengerId: Consider it the serial number of the data records. Each passenger is given this unique id in the data frame.
#Survived: 0 means not survived, 1 means the passenger survived
#Pclass: Passenger Class (1-1st class, 2-2nd Class, 3- 3rd Class)
#Sex: Gender of the passenger
#Age: Age of the passenger
#Fare: Ticket price that the passenger paid


#Data Frame Subsetting
#Selecting only PassengerId and Name columns
select(titanic, c(PassengerId, Name))
#Selecting only those records where sex is male
filter(titanic, Sex == "male")
#Selecting the record of only those males who survived
filter(titanic, Sex == "male" & Survived == 1)

#Applying functions over rows and columns
# calculating sum of Age and Fare column (Column wise)
apply(select(titanic, Age, Fare), 2, sum, na.rm = TRUE)
# calculating sum of Age and Fare column (Row wise)
apply(select(titanic, Age, Fare), 1, sum, na.rm = TRUE)

#Grouping and applying aggregate functions

titanic %>% group_by(Sex) %>% summarise(mean(Fare))

#Data Frame Pivoting

x <- select(titanic, Sex, Pclass, Survived)
y <- x %>% group_by(Sex, Pclass) %>% summarise(sum(Survived))
pt <- dcast(y, Sex~Pclass)
pt

#Data Frame Melting

melt(pt, variable.name = "Pclass")

#Merging multiple data frames together
#We will first create the other data frame that we will merge with the titanic data frame
y = data.frame(Pclass = c(1,2,4), Class = c("1st Class", "2nd Class", "4th Class"))
#Left Join:
merge(titanic, y, by = "Pclass", all.x = TRUE)
#Right Join:
merge(titanic, y, by = "Pclass", all.y = TRUE)
#Inner Join:
merge(titanic, y, by = "Pclass")
#Outer/Full Join:
merge(titanic, y, by = "Pclass", all.x = TRUE, all.y = TRUE)

