class Person:
  def __init__(self, name, age):
    self.name = name
    self.age = age

  def myfunc(self):
    print("Hello my name is " + self.name)

def runperson():
	p1 = Person("John", 36)
	p1.myfunc()
	
if __name__ == "__main__":
    runperson()
    
# python3 base_transformer_shanker/person.py
