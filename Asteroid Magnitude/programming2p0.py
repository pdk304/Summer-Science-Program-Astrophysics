import numpy as np

fruits = np.array([["Apple","Banana","Blueberry","Cherry"],
["Coconut","Grapefruit","Kumquat","Mango"],
["Nectarine","Orange","Tangerine","Pomegranate"],
["Lemon","Raspberry","Strawberry","Tomato"]])

#part(a)
print("part(a)")
print(fruits[3,3])

#part(b)
print("part(b)")
print(fruits[1:3,1:3])

#part(c)
print("part(c)")
print(fruits[0:3:2,0:4])

#part(d)
print("part(d)")
print(np.flip(np.flip(fruits[1:3,1:3],0),1))

#part(e)
print("part(e)")
copy = fruits.copy()
copy[:,0:1] = fruits[:,3:4]
copy[:,3:4] = fruits[:,0:1]
print(copy)

fruits[:,:] = "Sliced!"
print(fruits)
